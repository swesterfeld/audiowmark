/*
 * Copyright (C) 2018-2020 Stefan Westerfeld
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include <string.h>
#include <stdio.h>

#include <regex>

#include "utils.hh"
#include "mpegts.hh"
#include "wavdata.hh"
#include "wmcommon.hh"
#include "sfinputstream.hh"

extern "C" {
#include <libavformat/avformat.h>
#include <libavutil/opt.h>
#include <libswresample/swresample.h>
#include <libavutil/avassert.h>
#include <libavutil/timestamp.h>
#undef av_err2str
#define av_err2str(errnum) av_make_error_string((char*)__builtin_alloca(AV_ERROR_MAX_STRING_SIZE), AV_ERROR_MAX_STRING_SIZE, errnum)
}

using std::string;
using std::regex;
using std::vector;
using std::map;
using std::min;

Error
xsystem (const string& cmd)
{
  info ("+++ %s\n", cmd.c_str());
  int rc = system (cmd.c_str());
  int exit_status = WEXITSTATUS (rc);
  if (exit_status != 0)
    {
      error ("audiowmark: failed to execute command:\n%s\n", cmd.c_str());
      return Error (string_printf ("system failed / exit status %d", exit_status));
    }
  return Error::Code::NONE;
}

Error
ff_decode (const string& filename, WavData& out_wav_data)
{
  FILE *tmp_file = tmpfile();
  ScopedFile tmp_file_s (tmp_file);
  string tmp_file_name = string_printf ("/dev/fd/%d", fileno (tmp_file));

  FILE *input_tmp_file = tmpfile();
  ScopedFile input_tmp_file_s (input_tmp_file);
  string input_tmp_file_name = string_printf ("/dev/fd/%d", fileno (input_tmp_file));

  // write current ts
  FILE *main = fopen (filename.c_str(), "r");
  ScopedFile main_s (main);
  int c;
  while ((c = fgetc (main)) >= 0)
    fputc (c, input_tmp_file);

  fflush (input_tmp_file);
  string cmd = string_printf ("ffmpeg -v error -y -f mpegts -i %s -f wav %s", input_tmp_file_name.c_str(), tmp_file_name.c_str());
  Error err = xsystem (cmd.c_str());
  if (err)
    return err;

  err = out_wav_data.load (tmp_file_name);
  return err;
}

/*------------------------------- start code from ffmpeg...muxing.c -------------------------------*/

/*
 * Copyright (c) 2003 Fabrice Bellard
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
 * THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE.
 */

// a wrapper around a single output AVStream
typedef struct OutputStream {
    AVStream *st;
    AVCodecContext *enc;

    /* pts of the next frame that will be generated */
    int64_t next_pts;
    int samples_count;

    AVFrame *frame;
    AVFrame *tmp_frame;

    float t, tincr, tincr2;

    struct SwsContext *sws_ctx;
    struct SwrContext *swr_ctx;
} OutputStream;


/* Add an output stream. */
static void add_stream(OutputStream *ost, AVFormatContext *oc,
                       AVCodec **codec,
                       enum AVCodecID codec_id)
{
    AVCodecContext *c;
    int i;

    /* find the encoder */
    *codec = avcodec_find_encoder(codec_id);
    if (!(*codec)) {
        fprintf(stderr, "Could not find encoder for '%s'\n",
                avcodec_get_name(codec_id));
        exit(1);
    }

    ost->st = avformat_new_stream(oc, NULL);
    if (!ost->st) {
        fprintf(stderr, "Could not allocate stream\n");
        exit(1);
    }
    ost->st->id = oc->nb_streams-1;
    c = avcodec_alloc_context3(*codec);
    if (!c) {
        fprintf(stderr, "Could not alloc an encoding context\n");
        exit(1);
    }
    ost->enc = c;

    switch ((*codec)->type) {
    case AVMEDIA_TYPE_AUDIO:
        c->sample_fmt  = (*codec)->sample_fmts ?
            (*codec)->sample_fmts[0] : AV_SAMPLE_FMT_FLTP;
        c->bit_rate    = 64000;
        c->sample_rate = 44100;
        if ((*codec)->supported_samplerates) {
            c->sample_rate = (*codec)->supported_samplerates[0];
            for (i = 0; (*codec)->supported_samplerates[i]; i++) {
                if ((*codec)->supported_samplerates[i] == 44100)
                    c->sample_rate = 44100;
            }
        }
        c->channels        = av_get_channel_layout_nb_channels(c->channel_layout);
        c->channel_layout = AV_CH_LAYOUT_STEREO;
        if ((*codec)->channel_layouts) {
            c->channel_layout = (*codec)->channel_layouts[0];
            for (i = 0; (*codec)->channel_layouts[i]; i++) {
                if ((*codec)->channel_layouts[i] == AV_CH_LAYOUT_STEREO)
                    c->channel_layout = AV_CH_LAYOUT_STEREO;
            }
        }
        c->channels        = av_get_channel_layout_nb_channels(c->channel_layout);
        ost->st->time_base = (AVRational){ 1, c->sample_rate };
        break;

#if 0
    case AVMEDIA_TYPE_VIDEO:
        c->codec_id = codec_id;

        c->bit_rate = 400000;
        /* Resolution must be a multiple of two. */
        c->width    = 352;
        c->height   = 288;
        /* timebase: This is the fundamental unit of time (in seconds) in terms
         * of which frame timestamps are represented. For fixed-fps content,
         * timebase should be 1/framerate and timestamp increments should be
         * identical to 1. */
        ost->st->time_base = (AVRational){ 1, STREAM_FRAME_RATE };
        c->time_base       = ost->st->time_base;

        c->gop_size      = 12; /* emit one intra frame every twelve frames at most */
        c->pix_fmt       = STREAM_PIX_FMT;
        if (c->codec_id == AV_CODEC_ID_MPEG2VIDEO) {
            /* just for testing, we also add B-frames */
            c->max_b_frames = 2;
        }
        if (c->codec_id == AV_CODEC_ID_MPEG1VIDEO) {
            /* Needed to avoid using macroblocks in which some coeffs overflow.
             * This does not happen with normal video, it just happens here as
             * the motion of the chroma plane does not match the luma plane. */
            c->mb_decision = 2;
        }
    break;
#endif
    default:
        break;
    }

    /* Some formats want stream headers to be separate. */
    if (oc->oformat->flags & AVFMT_GLOBALHEADER)
        c->flags |= AV_CODEC_FLAG_GLOBAL_HEADER;
}


static AVFrame *alloc_audio_frame(enum AVSampleFormat sample_fmt,
                                  uint64_t channel_layout,
                                  int sample_rate, int nb_samples)
{
    AVFrame *frame = av_frame_alloc();
    int ret;

    if (!frame) {
        fprintf(stderr, "Error allocating an audio frame\n");
        exit(1);
    }

    frame->format = sample_fmt;
    frame->channel_layout = channel_layout;
    frame->sample_rate = sample_rate;
    frame->nb_samples = nb_samples;

    if (nb_samples) {
        ret = av_frame_get_buffer(frame, 0);
        if (ret < 0) {
            fprintf(stderr, "Error allocating an audio buffer\n");
            exit(1);
        }
    }

    return frame;
}


static void open_audio(AVFormatContext *oc, AVCodec *codec, OutputStream *ost, AVDictionary *opt_arg)
{
    AVCodecContext *c;
    int nb_samples;
    int ret;
    AVDictionary *opt = NULL;

    c = ost->enc;

    /* open it */
    av_dict_copy(&opt, opt_arg, 0);
    ret = avcodec_open2(c, codec, &opt);
    av_dict_free(&opt);
    if (ret < 0) {
        fprintf(stderr, "Could not open audio codec: %s\n", av_err2str(ret));
        exit(1);
    }

    /* init signal generator */
    ost->t     = 0;
    ost->tincr = 2 * M_PI * 110.0 / c->sample_rate;
    /* increment frequency by 110 Hz per second */
    ost->tincr2 = 2 * M_PI * 110.0 / c->sample_rate / c->sample_rate;

    if (c->codec->capabilities & AV_CODEC_CAP_VARIABLE_FRAME_SIZE)
        nb_samples = 10000;
    else
        nb_samples = c->frame_size;

    ost->frame     = alloc_audio_frame(c->sample_fmt, c->channel_layout,
                                       c->sample_rate, nb_samples);
    ost->tmp_frame = alloc_audio_frame(AV_SAMPLE_FMT_S16, c->channel_layout,
                                       c->sample_rate, nb_samples);

    /* copy the stream parameters to the muxer */
    ret = avcodec_parameters_from_context(ost->st->codecpar, c);
    if (ret < 0) {
        fprintf(stderr, "Could not copy the stream parameters\n");
        exit(1);
    }

    /* create resampler context */
        ost->swr_ctx = swr_alloc();
        if (!ost->swr_ctx) {
            fprintf(stderr, "Could not allocate resampler context\n");
            exit(1);
        }

        /* set options */
        av_opt_set_int       (ost->swr_ctx, "in_channel_count",   c->channels,       0);
        av_opt_set_int       (ost->swr_ctx, "in_sample_rate",     c->sample_rate,    0);
        av_opt_set_sample_fmt(ost->swr_ctx, "in_sample_fmt",      AV_SAMPLE_FMT_S16, 0);
        av_opt_set_int       (ost->swr_ctx, "out_channel_count",  c->channels,       0);
        av_opt_set_int       (ost->swr_ctx, "out_sample_rate",    c->sample_rate,    0);
        av_opt_set_sample_fmt(ost->swr_ctx, "out_sample_fmt",     c->sample_fmt,     0);

        /* initialize the resampling context */
        if ((ret = swr_init(ost->swr_ctx)) < 0) {
            fprintf(stderr, "Failed to initialize the resampling context\n");
            exit(1);
        }
}

#define STREAM_DURATION 0.1

/* Prepare a 16 bit dummy audio frame of 'frame_size' samples and
 * 'nb_channels' channels. */
static AVFrame *get_audio_frame(OutputStream *ost)
{
    AVFrame *frame = ost->tmp_frame;
    int j, i, v;
    int16_t *q = (int16_t*)frame->data[0];

    /* check if we want to generate more frames */
    if (av_compare_ts(ost->next_pts, ost->enc->time_base,
                      STREAM_DURATION, (AVRational){ 1, 1 }) > 0)
        return NULL;

    for (j = 0; j <frame->nb_samples; j++) {
        v = (int)(sin(ost->t) * 10000);
        for (i = 0; i < ost->enc->channels; i++)
            *q++ = v;
        ost->t     += ost->tincr;
        ost->tincr += ost->tincr2;
    }

    frame->pts = ost->next_pts;
    ost->next_pts  += frame->nb_samples;

    return frame;
}


static void log_packet(const AVFormatContext *fmt_ctx, const AVPacket *pkt)
{
    AVRational *time_base = &fmt_ctx->streams[pkt->stream_index]->time_base;

#if 0 // FIXME
    printf("pts:%s pts_time:%s dts:%s dts_time:%s duration:%s duration_time:%s stream_index:%d\n",
           av_ts2str(pkt->pts), av_ts2timestr(pkt->pts, time_base),
           av_ts2str(pkt->dts), av_ts2timestr(pkt->dts, time_base),
           av_ts2str(pkt->duration), av_ts2timestr(pkt->duration, time_base),
           pkt->stream_index);
#endif
}


static int write_frame(AVFormatContext *fmt_ctx, const AVRational *time_base, AVStream *st, AVPacket *pkt)
{
    /* rescale output packet timestamp values from codec to stream timebase */
    av_packet_rescale_ts(pkt, *time_base, st->time_base);
    pkt->stream_index = st->index;

    /* Write the compressed frame to the media file. */
    log_packet(fmt_ctx, pkt);
    return av_interleaved_write_frame(fmt_ctx, pkt);
}


/*
 * encode one audio frame and send it to the muxer
 * return 1 when encoding is finished, 0 otherwise
 */
static int write_audio_frame(AVFormatContext *oc, OutputStream *ost)
{
    AVCodecContext *c;
    AVPacket pkt = { 0 }; // data and size must be 0;
    AVFrame *frame;
    int ret;
    int got_packet;
    int dst_nb_samples;

    av_init_packet(&pkt);
    c = ost->enc;

    frame = get_audio_frame(ost);

    if (frame) {
        /* convert samples from native format to destination codec format, using the resampler */
            /* compute destination number of samples */
            dst_nb_samples = av_rescale_rnd(swr_get_delay(ost->swr_ctx, c->sample_rate) + frame->nb_samples,
                                            c->sample_rate, c->sample_rate, AV_ROUND_UP);
            av_assert0(dst_nb_samples == frame->nb_samples);

        /* when we pass a frame to the encoder, it may keep a reference to it
         * internally;
         * make sure we do not overwrite it here
         */
        ret = av_frame_make_writable(ost->frame);
        if (ret < 0)
            exit(1);

        /* convert to destination format */
        ret = swr_convert(ost->swr_ctx,
                          ost->frame->data, dst_nb_samples,
                          (const uint8_t **)frame->data, frame->nb_samples);
        if (ret < 0) {
            fprintf(stderr, "Error while converting\n");
            exit(1);
        }
        frame = ost->frame;

        frame->pts = av_rescale_q(ost->samples_count, (AVRational){1, c->sample_rate}, c->time_base);
        ost->samples_count += dst_nb_samples;
    }

    ret = avcodec_encode_audio2(c, &pkt, frame, &got_packet);
    if (ret < 0) {
        fprintf(stderr, "Error encoding audio frame: %s\n", av_err2str(ret));
        exit(1);
    }

    if (got_packet) {
        ret = write_frame(oc, &c->time_base, ost->st, &pkt);
        if (ret < 0) {
            fprintf(stderr, "Error while writing audio frame: %s\n",
                    av_err2str(ret));
            exit(1);
        }
    }

    return (frame || got_packet) ? 0 : 1;
}

static void close_stream(AVFormatContext *oc, OutputStream *ost)
{
    avcodec_free_context(&ost->enc);
    av_frame_free(&ost->frame);
    av_frame_free(&ost->tmp_frame);
    swr_free(&ost->swr_ctx);
}


Error
ff_encode (const WavData& wav_data, const string& filename, size_t start_pos, size_t cut_start, size_t cut_end, double pts_start)
{
  AVFormatContext *oc;
  avformat_alloc_output_context2 (&oc, NULL, "mpegts", NULL);
  if (!oc)
    return Error ("failed to alloc avformat output context");
  int ret = avio_open (&oc->pb, "pipe:1", AVIO_FLAG_WRITE);
  if (ret < 0)
    {
      error ("Could not open output: %s\n", av_err2str (ret));
      return Error ("open pipe failed");
    }

  OutputStream audio_st = { 0 };
  AVCodec *audio_codec;
  AVOutputFormat *fmt;
  AVDictionary *opt = nullptr;

  fmt = oc->oformat;
  add_stream (&audio_st, oc, &audio_codec, AV_CODEC_ID_AAC);
  open_audio(oc, audio_codec, &audio_st, opt);

  /* Write the stream header, if any. */
  ret = avformat_write_header (oc, &opt);
  if (ret < 0)
    {
      error ("Error occurred when writing output file: %s\n",  av_err2str(ret));
      return Error ("avformat_write_header failed\n");
    }
  av_dump_format(oc, 0, "pipe:1", 1);
  while (write_audio_frame(oc, &audio_st) == 0);
  av_write_trailer(oc);

  close_stream(oc, &audio_st);

  if (!(fmt->flags & AVFMT_NOFILE))
      /* Close the output file. */
      avio_closep(&oc->pb);

  /* free the stream */
  avformat_free_context(oc);

/*------------------------------- end code from ffmpeg...muxing.c -------------------------------*/


#if 0
  FILE *tmp_file = tmpfile();
  ScopedFile tmp_file_s (tmp_file);
  string tmp_file_name = string_printf ("/dev/fd/%d", fileno (tmp_file));

  Error err = wav_data.save (tmp_file_name);

  string cmd = string_printf ("ffmpeg -v error -y -i %s -f mpegts -c:a aac '%s'", tmp_file_name.c_str(), filename.c_str());
  err = xsystem (cmd);
  if (err)
    return err;

  double length_s = double (wav_data.n_values()) / wav_data.n_channels() / wav_data.sample_rate();

  // cut_start_s is corrected down to avoid cutting one frame more or than intended
  double cut_start_s = cut_start / double (wav_data.sample_rate());
  double cut_end_s = cut_end / double (wav_data.sample_rate());
  cmd = string_printf ("ffmpeg -v error -y -i '%s' -ss %.6f -t %.6f -f mpegts -output_ts_offset %f -muxdelay 0 -muxpreload 0 -c copy '%s-tcpy'",
                       filename.c_str(), cut_start_s - 1. / wav_data.sample_rate(), length_s - (cut_start_s + cut_end_s), pts_start, filename.c_str());
  err = xsystem (cmd);
  if (err)
    return err;

  cmd = string_printf ("mv '%s-tcpy' '%s'", filename.c_str(), filename.c_str());
  xsystem (cmd);
  if (err)
    return err;
#endif
  return Error::Code::NONE;
}

int
hls_embed_context (const string& in_dir, const string& out_dir, const string& filename, const string& audio_master)
{
  string in_name = in_dir + "/" + filename;
  FILE *in_file = fopen (in_name.c_str(), "r");
  ScopedFile in_file_s (in_file);

  if (!in_file)
    {
      error ("audiowmark: error opening input playlist %s\n", in_name.c_str());
      return 1;
    }

  string out_name = out_dir + "/" + filename;
  FILE *out_file = fopen (out_name.c_str(), "w");
  ScopedFile out_file_s (out_file);

  if (!out_file)
    {
      error ("audiowmark: error opening output playlist %s\n", out_name.c_str());
      return 1;
    }

  WavData audio_master_data;
  Error err = audio_master_data.load (audio_master);
  if (err)
    {
      error ("audiowmark: failed to load audio master: %s\n", audio_master.c_str());
      return 1;
    }

  struct Segment
  {
    string              name;
    size_t              size;
    map<string, string> vars;
  };
  vector<Segment> segments;
  char buffer[1024];
  int line = 1;
  const regex blank_re (R"(\s*(#.*)?)");
  while (fgets (buffer, 1024, in_file))
    {
      /* kill newline chars at end */
      int last = strlen (buffer) - 1;
      while (last > 0 && (buffer[last] == '\n' || buffer[last] == '\r'))
        buffer[last--] = 0;

      string s = buffer;

      std::smatch match;
      if (regex_match (s, blank_re))
        {
          /* blank line or comment */
          fprintf (out_file, "%s\n", s.c_str());
        }
      else
        {
          fprintf (out_file, "%s\n", s.c_str());
          Segment segment;
          segment.name = s;
          segments.push_back (segment);
        }
      line++;
    }
  size_t start_pos = 0;
  for (auto& segment : segments)
    {
      WavData out;
      Error err = ff_decode (in_dir + "/" + segment.name, out);
      if (err)
        {
          error ("audiowmark: hls: ff_decode failed: %s\n", err.message());
          return 1;
        }
      printf ("%d %zd\n", out.sample_rate(), out.n_values() / out.n_channels());
      segment.size = out.n_values() / out.n_channels();

      /* obtain pts for first frame */
      string cmd = string_printf ("ffprobe -v 0 -show_entries packet=pts_time %s/%s -of compact=p=0:nk=1 | grep '^[0-9]'", in_dir.c_str(), segment.name.c_str());
      FILE *pts = popen (cmd.c_str(), "r");
      char buffer[1024];
      if (fgets (buffer, 1024, pts))
        {
          if (strlen (buffer) && buffer[strlen (buffer) - 1] == '\n')
            buffer[strlen (buffer) - 1] = 0;
          segment.vars["pts_start"] = buffer;
        }
      fclose (pts);

      segment.vars["start_pos"] = string_printf ("%zd", start_pos);
      segment.vars["size"] = string_printf ("%zd", segment.size);

      start_pos += segment.size;
    }

  /* fill out next/prev size fields */
  for (size_t i = 0; i < segments.size(); i++)
    {
      if (i > 0)
        segments[i].vars["prev_size"] = string_printf ("%zd", segments[i - 1].size);
      else
        segments[i].vars["prev_size"] = "0";

      if (i + 1 < segments.size())
        segments[i].vars["next_size"] = string_printf ("%zd", segments[i + 1].size);
      else
        segments[i].vars["next_size"] = "0";
    }

  /* write audio segments with context */
  for (auto& segment : segments)
    {
      /* write a part of audio master here */
      const size_t prev_size = atoi (segment.vars["prev_size"].c_str());
      const size_t next_size = atoi (segment.vars["next_size"].c_str());

      const size_t start_point = atoi (segment.vars["start_pos"].c_str()) - prev_size;
      const size_t end_point = start_point + prev_size + segment.size + next_size;

      vector<float> out_signal (audio_master_data.samples().begin() + start_point * audio_master_data.n_channels(),
                                audio_master_data.samples().begin() + end_point * audio_master_data.n_channels());
      WavData out_wav_data (out_signal, audio_master_data.n_channels(), audio_master_data.sample_rate(), audio_master_data.bit_depth());
      err = out_wav_data.save (out_dir + "/" + segment.name + ".wav");
    }

  for (size_t i = 0; i < segments.size(); i++)
    {
      TSWriter writer;

      writer.append_file ("full.wav", out_dir + "/" + segments[i].name + ".wav");
      writer.append_vars ("vars", segments[i].vars);
      writer.process (in_dir + "/" + segments[i].name, out_dir + "/" + segments[i].name);
    }
  return 0;
}

class WDInputStream : public AudioInputStream
{
  WavData *wav_data;
  size_t   read_pos = 0;
public:
  WDInputStream (WavData *wav_data) :
    wav_data (wav_data)
  {
  }
  int
  bit_depth() const override
  {
    return wav_data->bit_depth();
  }
  int
  sample_rate() const override
  {
    return wav_data->sample_rate();
  }
  int
  n_channels() const override
  {
    return wav_data->n_channels();
  }
  size_t
  n_frames() const override
  {
    return wav_data->n_values() / wav_data->n_channels();
  }
  Error
  read_frames (std::vector<float>& samples, size_t count) override
  {
    size_t read_count = min (n_frames() - read_pos, count);

    const auto& wsamples = wav_data->samples();
    samples.assign (wsamples.begin() + read_pos * n_channels(), wsamples.begin() + (read_pos + read_count) * n_channels());

    read_pos += read_count;

    return Error::Code::NONE;
  }
};

class WDOutputStream : public AudioOutputStream
{
  WavData *wav_data;
  vector<float> samples;
public:
  WDOutputStream (WavData *wav_data) :
    wav_data (wav_data)
  {
  }
  int
  bit_depth() const override
  {
    return wav_data->bit_depth();
  }
  int
  sample_rate() const override
  {
    return wav_data->sample_rate();
  }
  int
  n_channels() const override
  {
    return wav_data->n_channels();
  }
  Error
  write_frames (const std::vector<float>& frames)
  {
    samples.insert (samples.end(), frames.begin(), frames.end());
    return Error::Code::NONE;
  }
  Error
  close()
  {
    wav_data->set_samples (samples); // only do this once at end for performance reasons
    return Error::Code::NONE;
  }
};

int
mark_zexpand (WavData& wav_data, size_t zero_frames, const string& bits)
{
  WDInputStream in_stream (&wav_data);

  WavData wav_data_out ({ /* no samples */ }, wav_data.n_channels(), wav_data.sample_rate(), wav_data.bit_depth());
  WDOutputStream out_stream (&wav_data_out);

  int rc = add_stream_watermark (&in_stream, &out_stream, bits, zero_frames);
  if (rc != 0)
    return rc;

  wav_data.set_samples (wav_data_out.samples());

  return 0;
}

int
hls_mark (const string& infile, const string& outfile, const string& bits)
{
  double start_time = get_time();

  TSReader reader;

  Error err = reader.load (infile);
  if (err)
    {
      error ("hls_mark: %s\n", err.message());
      return 1;
    }
  printf ("hls_elapsed_load %f\n", (get_time() - start_time) * 1000 /* ms */);
  double start_time1 = get_time();

  const TSReader::Entry *full_wav = reader.find ("full.wav");
  if (!full_wav)
    {
      error ("hls_mark: no embedded context found in %s\n", infile.c_str());
      return 1;
    }

  SFInputStream in_stream;
  err = in_stream.open (&full_wav->data);
  if (err)
    {
      error ("hls_mark: %s\n", err.message());
      return 1;
    }

  for (auto entry : reader.entries())
    printf ("%s %zd\n", entry.filename.c_str(), entry.data.size());

  map<string, string> vars = reader.parse_vars ("vars");
  for (auto kv : vars)
    printf ("|| %s=%s\n", kv.first.c_str(), kv.second.c_str());

  size_t start_pos = atoi (vars["start_pos"].c_str());
  size_t prev_size = atoi (vars["prev_size"].c_str());
  size_t next_size = atoi (vars["next_size"].c_str());
  double pts_start = atof (vars["pts_start"].c_str());
  size_t next_ctx = min<size_t> (1024 * 3, next_size);
  size_t prev_ctx = min<size_t> (1024 * 3, prev_size);

  printf ("hls_time_elapsed_decode %f\n", (get_time() - start_time1) * 1000 /* ms */);
  start_time1 = get_time();

  WavData wav_data ({ /* no samples */ }, in_stream.n_channels(), in_stream.sample_rate(), in_stream.bit_depth());
  WDOutputStream out_stream (&wav_data);

  int zrc = add_stream_watermark (&in_stream, &out_stream, bits, start_pos - prev_size);
  if (zrc != 0)
    return zrc;

  /* erase extra samples caused by concatting with prev.ts */
  auto samples = wav_data.samples();
  samples.erase (samples.begin(), samples.begin() + (prev_size - prev_ctx) * wav_data.n_channels());
  samples.erase (samples.end() - (next_size - next_ctx) * wav_data.n_channels(), samples.end());
  wav_data.set_samples (samples);

  printf ("hls_time_elapsed_mark %f\n", (get_time() - start_time1) * 1000 /* ms */);
  start_time1 = get_time();

  err = ff_encode (wav_data, outfile, start_pos, start_pos == 0 ? 1024 : prev_ctx, next_ctx, pts_start);
  if (err)
    {
      error ("hls_mark: %s\n", err.message());
      return 1;
    }

  printf ("hls_time_elapsed_aac_enc %f\n", (get_time() - start_time1) * 1000 /* ms */);

  double end_time = get_time();
  printf ("hls_time %f %f\n", start_pos / double (wav_data.sample_rate()), (end_time - start_time) * 1000 /* ms */);

  return 0;
}

int
test_seek (const string& in, const string& out, int pos, const string& bits)
{
  vector<float> samples;
  WavData wav_data;
  Error err = wav_data.load (in);
  if (err)
    {
      error ("load error: %s\n", err.message());
      return 1;
    }

  samples = wav_data.samples();
  samples.erase (samples.begin(), samples.begin() + pos * wav_data.n_channels());
  wav_data.set_samples (samples);

  int rc = mark_zexpand (wav_data, pos, bits);
  if (rc != 0)
    {
      return rc;
    }

  samples = wav_data.samples();
  samples.insert (samples.begin(), pos * wav_data.n_channels(), 0);
  wav_data.set_samples (samples);

  err = wav_data.save (out);
  if (err)
    {
      error ("save error: %s\n", err.message());
      return 1;
    }
  return 0;
}

int
seek_perf (int sample_rate, double seconds)
{
  vector<float> samples (100);
  WavData wav_data (samples, 2, sample_rate, 16);

  double start_time = get_time();

  int rc = mark_zexpand (wav_data, seconds * sample_rate, "0c");
  if (rc != 0)
    return rc;

  double end_time = get_time();

  printf ("\n\n");
  printf ("total time %7.3f sec\n", end_time - start_time);
  printf ("per second %7.3f ms\n", (end_time - start_time) / seconds * 1000);

  return 0;
}

int
main (int argc, char **argv)
{
  if (argc == 6 && strcmp (argv[1], "hls-embed-context") == 0)
    {
      printf ("hls-embed-context: in_dir=%s out_dir=%s m3u8=%s audio_master=%s\n", argv[2], argv[3], argv[4], argv[5]);
      return hls_embed_context (argv[2], argv[3], argv[4], argv[5]);
    }
  else if (argc == 5 && strcmp (argv[1], "hls-mark") == 0)
    {
      return hls_mark (argv[2], argv[3], argv[4]);
    }
  else if (argc == 6 && strcmp (argv[1], "test-seek") == 0)
    {
      return test_seek (argv[2], argv[3], atoi (argv[4]), argv[5]);
    }
  else if (argc == 4 && strcmp (argv[1], "seek-perf") == 0)
    {
      return seek_perf (atoi (argv[2]), atof (argv[3]));
    }
  else if (argc == 4 && strcmp (argv[1], "ff-decode") == 0)
    {
      WavData wd;
      Error err = ff_decode (argv[2], wd);
      if (err)
        {
          error ("audiowmark: hls: ff_decode failed: %s\n", err.message());
          return 1;
        }
      err = wd.save (argv[3]);
      if (err)
        {
          error ("audiowmark: hls: save failed: %s\n", err.message());
          return 1;
        }
      return 0;
    }
  else
    {
      error ("testhls: error parsing command line arguments\n");
      return 1;
    }
}

