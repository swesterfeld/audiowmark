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

#include "hlsoutputstream.hh"

/* HLSOutputStream is based on code from ffmpeg: doc/examples/muxing.c */

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

using std::vector;
using std::string;
using std::min;

HLSOutputStream::HLSOutputStream (int n_channels, int sample_rate, int bit_depth) :
  m_bit_depth (bit_depth),
  m_sample_rate (sample_rate),
  m_n_channels (n_channels),
  m_audio_buffer (n_channels)
{
}

/* Add an output stream. */
void
HLSOutputStream::add_stream (AVCodec **codec, enum AVCodecID codec_id)
{
  /* find the encoder */
  *codec = avcodec_find_encoder (codec_id);
  if (!(*codec))
    {
      fprintf(stderr, "Could not find encoder for '%s'\n", avcodec_get_name (codec_id));
      exit(1);
    }

  m_st = avformat_new_stream (m_fmt_ctx, NULL);
  if (!m_st)
    {
      fprintf (stderr, "Could not allocate stream\n");
      exit(1);
    }
  m_st->id = m_fmt_ctx->nb_streams - 1;

  m_enc = avcodec_alloc_context3 (*codec);
  if (!m_enc)
    {
      fprintf (stderr, "Could not alloc an encoding context\n");
      exit(1);
    }

  if ((*codec)->type != AVMEDIA_TYPE_AUDIO)
    {
      error ("HLSOutputStream: codec type must be audio");
      exit (1);
    }

  m_enc->sample_fmt  = (*codec)->sample_fmts ? (*codec)->sample_fmts[0] : AV_SAMPLE_FMT_FLTP;
  m_enc->bit_rate    = 128000;
  m_enc->sample_rate = 44100;
  if ((*codec)->supported_samplerates)
    {
      m_enc->sample_rate = (*codec)->supported_samplerates[0];
        for (int i = 0; (*codec)->supported_samplerates[i]; i++)
          {
            if ((*codec)->supported_samplerates[i] == 44100)
              m_enc->sample_rate = 44100;
          }
    }
  m_enc->channels       = av_get_channel_layout_nb_channels (m_enc->channel_layout);
  m_enc->channel_layout = AV_CH_LAYOUT_STEREO;
  if ((*codec)->channel_layouts)
    {
      m_enc->channel_layout = (*codec)->channel_layouts[0];
      for (int i = 0; (*codec)->channel_layouts[i]; i++)
        {
          if ((*codec)->channel_layouts[i] == AV_CH_LAYOUT_STEREO)
              m_enc->channel_layout = AV_CH_LAYOUT_STEREO;
        }
    }
  m_enc->channels     = av_get_channel_layout_nb_channels (m_enc->channel_layout);
  m_st->time_base = (AVRational){ 1, m_enc->sample_rate };

  /* Some formats want stream headers to be separate. */
  if (m_fmt_ctx->oformat->flags & AVFMT_GLOBALHEADER)
    m_enc->flags |= AV_CODEC_FLAG_GLOBAL_HEADER;
}


AVFrame *
HLSOutputStream::alloc_audio_frame (AVSampleFormat sample_fmt, uint64_t channel_layout, int sample_rate, int nb_samples)
{
  AVFrame *frame = av_frame_alloc();

  if (!frame)
    {
      fprintf (stderr, "Error allocating an audio frame\n");
      exit(1);
    }

  frame->format = sample_fmt;
  frame->channel_layout = channel_layout;
  frame->sample_rate = sample_rate;
  frame->nb_samples = nb_samples;

  if (nb_samples)
    {
      int ret = av_frame_get_buffer (frame, 0);
      if (ret < 0)
        {
          fprintf (stderr, "Error allocating an audio buffer\n");
          exit(1);
        }
    }

  return frame;
}


void
HLSOutputStream::open_audio (AVCodec *codec, AVDictionary *opt_arg)
{
  int nb_samples;
  int ret;
  AVDictionary *opt = NULL;

  /* open it */
  av_dict_copy (&opt, opt_arg, 0);
  ret = avcodec_open2 (m_enc, codec, &opt);
  av_dict_free (&opt);
  if (ret < 0)
    {
      fprintf(stderr, "Could not open audio codec: %s\n", av_err2str(ret));
      exit(1);
    }

  if (m_enc->codec->capabilities & AV_CODEC_CAP_VARIABLE_FRAME_SIZE)
    nb_samples = 10000;
  else
    nb_samples = m_enc->frame_size;

  m_frame     = alloc_audio_frame (m_enc->sample_fmt, m_enc->channel_layout, m_enc->sample_rate, nb_samples);
  m_tmp_frame = alloc_audio_frame (AV_SAMPLE_FMT_S16, m_enc->channel_layout, m_enc->sample_rate, nb_samples);

  /* copy the stream parameters to the muxer */
  ret = avcodec_parameters_from_context (m_st->codecpar, m_enc);
  if (ret < 0)
    {
      fprintf(stderr, "Could not copy the stream parameters\n");
      exit(1);
    }

  /* create resampler context */
  m_swr_ctx = swr_alloc();
  if (!m_swr_ctx)
    {
      fprintf(stderr, "Could not allocate resampler context\n");
      exit(1);
    }

  /* set options */
  av_opt_set_int        (m_swr_ctx, "in_channel_count",   m_enc->channels,       0);
  av_opt_set_int        (m_swr_ctx, "in_sample_rate",     m_enc->sample_rate,    0);
  av_opt_set_sample_fmt (m_swr_ctx, "in_sample_fmt",      AV_SAMPLE_FMT_S16,     0);
  av_opt_set_int        (m_swr_ctx, "out_channel_count",  m_enc->channels,       0);
  av_opt_set_int        (m_swr_ctx, "out_sample_rate",    m_enc->sample_rate,    0);
  av_opt_set_sample_fmt (m_swr_ctx, "out_sample_fmt",     m_enc->sample_fmt,     0);

  /* initialize the resampling context */
  if ((ret = swr_init(m_swr_ctx)) < 0)
    {
      fprintf(stderr, "Failed to initialize the resampling context\n");
      exit(1);
    }
}

/* Prepare a 16 bit dummy audio frame of 'frame_size' samples and
 * 'nb_channels' channels. */
AVFrame *
HLSOutputStream::get_audio_frame()
{
  AVFrame *frame = m_tmp_frame;
  int16_t *q = (int16_t*)frame->data[0];

  if (m_audio_buffer.can_read_frames() < size_t (frame->nb_samples))
    return NULL;

  vector<float> samples = m_audio_buffer.read_frames (frame->nb_samples);

  size_t t = 0;
  for (int j = 0; j < frame->nb_samples; j++)
    {
      for (int i = 0; i < m_enc->channels; i++)
        {
          if (t < samples.size())
            {
              *q++ = (int)(samples[t] * 32768);
              t++;
            }
          else
            *q++ = 0;
        }
    }

  frame->pts = m_next_pts;
  m_next_pts  += frame->nb_samples;

  return frame;
}


int
HLSOutputStream::write_frame (const AVRational *time_base, AVStream *st, AVPacket *pkt)
{
  /* rescale output packet timestamp values from codec to stream timebase */
  av_packet_rescale_ts (pkt, *time_base, st->time_base);
  pkt->stream_index = st->index;

  /* Write the compressed frame to the media file. */
  return av_interleaved_write_frame (m_fmt_ctx, pkt);
}


/*
 * encode one audio frame and send it to the muxer
 * return 1 when encoding is finished, 0 otherwise
 */
int
HLSOutputStream::write_audio_frame()
{
  AVPacket pkt = { 0 }; // data and size must be 0;
  AVFrame *frame;
  int ret;
  int got_packet;

  av_init_packet (&pkt);

  frame = get_audio_frame();
  if (frame)
    {
      /* convert samples from native format to destination codec format, using the resampler */

      /* compute destination number of samples */
      int dst_nb_samples = av_rescale_rnd (swr_get_delay (m_swr_ctx, m_enc->sample_rate) + frame->nb_samples,
                                           m_enc->sample_rate, m_enc->sample_rate, AV_ROUND_UP);
      av_assert0 (dst_nb_samples == frame->nb_samples);

      /* when we pass a frame to the encoder, it may keep a reference to it
       * internally;
       * make sure we do not overwrite it here
       */
      ret = av_frame_make_writable (m_frame);
      if (ret < 0)
        exit(1);

      /* convert to destination format */
      ret = swr_convert (m_swr_ctx,
                         m_frame->data, dst_nb_samples,
                         (const uint8_t **)frame->data, frame->nb_samples);
      if (ret < 0)
        {
          fprintf (stderr, "Error while converting\n");
          exit(1);
        }
      frame = m_frame;

      frame->pts = av_rescale_q (m_samples_count + m_start_pos, (AVRational){1, m_enc->sample_rate}, m_enc->time_base);
      m_samples_count += dst_nb_samples;
    }

  ret = avcodec_encode_audio2 (m_enc, &pkt, frame, &got_packet);
  if (ret < 0)
    {
      fprintf(stderr, "Error encoding audio frame: %s\n", av_err2str(ret));
      exit(1);
    }

  if (got_packet)
    {
      if (m_cut_aac_frames)
        {
          m_cut_aac_frames--;
        }
      else if (m_keep_aac_frames)
        {
          ret = write_frame (&m_enc->time_base, m_st, &pkt);
          if (ret < 0)
            {
              fprintf(stderr, "Error while writing audio frame: %s\n",
                      av_err2str(ret));
              exit(1);
            }
          m_keep_aac_frames--;
        }
    }

  return (frame || got_packet) ? 0 : 1;
}

void
HLSOutputStream::close_stream()
{
  avcodec_free_context (&m_enc);
  av_frame_free (&m_frame);
  av_frame_free (&m_tmp_frame);
  swr_free (&m_swr_ctx);
}

Error
HLSOutputStream::open (const string& out_filename, size_t cut_aac_frames, size_t keep_aac_frames, double pts_start, size_t delete_input_start)
{
  avformat_alloc_output_context2 (&m_fmt_ctx, NULL, "mpegts", NULL);
  if (!m_fmt_ctx)
    return Error ("failed to alloc avformat output context");

  string filename = out_filename;
  if (filename == "-")
    filename = "pipe:1";

  int ret = avio_open (&m_fmt_ctx->pb, filename.c_str(), AVIO_FLAG_WRITE);
  if (ret < 0)
    {
      error ("Could not open output: %s\n", av_err2str (ret));
      return Error ("open hls output failed");
    }

  AVDictionary *opt = nullptr;
  AVCodec *audio_codec;
  add_stream (&audio_codec, AV_CODEC_ID_AAC);
  open_audio (audio_codec, opt);

  /* Write the stream header, if any. */
  ret = avformat_write_header (m_fmt_ctx, &opt);
  if (ret < 0)
    {
      error ("Error occurred when writing output file: %s\n",  av_err2str(ret));
      return Error ("avformat_write_header failed\n");
    }
  av_dump_format (m_fmt_ctx, 0, filename.c_str(), 1);

  m_delete_input_start = delete_input_start;
  m_cut_aac_frames = cut_aac_frames;
  m_keep_aac_frames = keep_aac_frames;

  // FIXME: correct?
  m_start_pos = pts_start * m_sample_rate - cut_aac_frames * 1024;
  m_start_pos += 1024;

  return Error::Code::NONE;
}

Error
HLSOutputStream::close()
{
  write(); // drain

  av_write_trailer (m_fmt_ctx);

  close_stream();

  /* Close the output file. */
  if (!(m_fmt_ctx->oformat->flags & AVFMT_NOFILE))
    avio_closep (&m_fmt_ctx->pb);

  /* free the stream */
  avformat_free_context (m_fmt_ctx);

  return Error::Code::NONE;
}

void
HLSOutputStream::write()
{
  while (write_audio_frame() == 0);
}

Error
HLSOutputStream::write_frames (const std::vector<float>& frames)
{
  // if we don't need any more aac frames, just throw away samples (save cpu cycles)
  if (m_keep_aac_frames == 0)
    return Error::Code::NONE;

  m_audio_buffer.write_frames (frames);

  size_t delete_input = min (m_delete_input_start, m_audio_buffer.can_read_frames());
  if (delete_input)
    {
      m_audio_buffer.read_frames (delete_input);
      m_delete_input_start -= delete_input;
    }

  while (m_audio_buffer.can_read_frames() >= 1024)
    {
      write_audio_frame();
    }
  return Error::Code::NONE;
}

int
HLSOutputStream::bit_depth() const
{
  return m_bit_depth;
}

int
HLSOutputStream::sample_rate() const
{
  return m_sample_rate;
}

int
HLSOutputStream::n_channels() const
{
  return m_n_channels;
}
