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

#undef av_err2str
#define av_err2str(errnum) av_make_error_string((char*)__builtin_alloca(AV_ERROR_MAX_STRING_SIZE), AV_ERROR_MAX_STRING_SIZE, errnum)

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
  av_log_set_level (AV_LOG_ERROR);
}

void
HLSOutputStream::set_bit_rate (int bit_rate)
{
  m_bit_rate = bit_rate;
}

void
HLSOutputStream::set_channel_layout (const string& channel_layout)
{
  m_channel_layout = channel_layout;
}

HLSOutputStream::~HLSOutputStream()
{
  close();
}

/* Add an output stream. */
Error
HLSOutputStream::add_stream (const AVCodec **codec, enum AVCodecID codec_id)
{
  /* find the encoder */
  *codec = avcodec_find_encoder (codec_id);
  if (!(*codec))
    return Error (string_printf ("could not find encoder for '%s'", avcodec_get_name (codec_id)));

  m_st = avformat_new_stream (m_fmt_ctx, NULL);
  if (!m_st)
    return Error ("could not allocate stream");

  m_st->id = m_fmt_ctx->nb_streams - 1;

  m_enc = avcodec_alloc_context3 (*codec);
  if (!m_enc)
    return Error ("could not alloc an encoding context");

  if ((*codec)->type != AVMEDIA_TYPE_AUDIO)
    return Error ("codec type must be audio");

  m_enc->sample_fmt  = (*codec)->sample_fmts ? (*codec)->sample_fmts[0] : AV_SAMPLE_FMT_FLTP;
  m_enc->bit_rate    = m_bit_rate;
  m_enc->sample_rate = m_sample_rate;
  if ((*codec)->supported_samplerates)
    {
      bool match = false;
      for (int i = 0; (*codec)->supported_samplerates[i]; i++)
        {
          if ((*codec)->supported_samplerates[i] == m_sample_rate)
            {
              m_enc->sample_rate = m_sample_rate;
              match = true;
            }
        }
      if (!match)
        return Error (string_printf ("no codec support for sample rate %d", m_sample_rate));
    }
  uint64_t want_layout = av_get_channel_layout (m_channel_layout.c_str());
  if (!want_layout)
    return Error (string_printf ("bad channel layout '%s'", m_channel_layout.c_str()));
  m_enc->channel_layout = want_layout;
  if ((*codec)->channel_layouts)
    {
      m_enc->channel_layout = (*codec)->channel_layouts[0];
      for (int i = 0; (*codec)->channel_layouts[i]; i++)
        {
          if ((*codec)->channel_layouts[i] == want_layout)
              m_enc->channel_layout = want_layout;
        }
    }
  if (want_layout != m_enc->channel_layout)
    return Error (string_printf ("codec: unsupported channel layout '%s'", m_channel_layout.c_str()));
  m_enc->channels = av_get_channel_layout_nb_channels (m_enc->channel_layout);
  m_st->time_base = (AVRational){ 1, m_enc->sample_rate };

  /* Some formats want stream headers to be separate. */
  if (m_fmt_ctx->oformat->flags & AVFMT_GLOBALHEADER)
    m_enc->flags |= AV_CODEC_FLAG_GLOBAL_HEADER;

  return Error::Code::NONE;
}


AVFrame *
HLSOutputStream::alloc_audio_frame (AVSampleFormat sample_fmt, uint64_t channel_layout, int sample_rate, int nb_samples, Error& err)
{
  AVFrame *frame = av_frame_alloc();

  if (!frame)
    {
      err = Error ("error allocating an audio frame");
      return nullptr;
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
          err = Error ("Error allocating an audio buffer");
          return nullptr;
        }
    }

  return frame;
}


Error
HLSOutputStream::open_audio (const AVCodec *codec, AVDictionary *opt_arg)
{
  int nb_samples;
  int ret;
  AVDictionary *opt = NULL;

  /* open it */
  av_dict_copy (&opt, opt_arg, 0);
  ret = avcodec_open2 (m_enc, codec, &opt);
  av_dict_free (&opt);
  if (ret < 0)
    return Error (string_printf ("could not open audio codec: %s", av_err2str (ret)));

  if (m_enc->codec->capabilities & AV_CODEC_CAP_VARIABLE_FRAME_SIZE)
    nb_samples = 10000;
  else
    nb_samples = m_enc->frame_size;

  Error err;
  m_frame     = alloc_audio_frame (m_enc->sample_fmt, m_enc->channel_layout, m_enc->sample_rate, nb_samples, err);
  if (err)
    return err;

  m_tmp_frame = alloc_audio_frame (AV_SAMPLE_FMT_FLT, m_enc->channel_layout, m_enc->sample_rate, nb_samples, err);
  if (err)
    return err;

  m_tmp_pkt = av_packet_alloc();
  if (!m_tmp_pkt)
    return Error ("could not allocate AVPacket");

  /* copy the stream parameters to the muxer */
  ret = avcodec_parameters_from_context (m_st->codecpar, m_enc);
  if (ret < 0)
    return Error ("could not copy the stream parameters");

  /* create resampler context */
  m_swr_ctx = swr_alloc();
  if (!m_swr_ctx)
    return Error ("could not allocate resampler context");

  /* set options */
  av_opt_set_int        (m_swr_ctx, "in_channel_count",   m_enc->channels,       0);
  av_opt_set_int        (m_swr_ctx, "in_sample_rate",     m_enc->sample_rate,    0);
  av_opt_set_sample_fmt (m_swr_ctx, "in_sample_fmt",      AV_SAMPLE_FMT_FLT,     0);
  av_opt_set_int        (m_swr_ctx, "out_channel_count",  m_enc->channels,       0);
  av_opt_set_int        (m_swr_ctx, "out_sample_rate",    m_enc->sample_rate,    0);
  av_opt_set_sample_fmt (m_swr_ctx, "out_sample_fmt",     m_enc->sample_fmt,     0);

  /* initialize the resampling context */
  if ((ret = swr_init(m_swr_ctx)) < 0)
    return Error ("failed to initialize the resampling context");

  return Error::Code::NONE;
}

/* fill audio frame with samples from AudioBuffer */
AVFrame *
HLSOutputStream::get_audio_frame()
{
  AVFrame *frame = m_tmp_frame;

  if (m_audio_buffer.can_read_frames() < size_t (frame->nb_samples))
    return nullptr;

  vector<float> samples = m_audio_buffer.read_frames (frame->nb_samples);

  std::copy (samples.begin(), samples.end(), (float *)frame->data[0]);

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
 *   returns EncResult: OK, ERROR, DONE
 */
HLSOutputStream::EncResult
HLSOutputStream::write_audio_frame (Error& err)
{
  AVFrame *frame;
  int ret;

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
        {
          err = Error ("error making frame writable");
          return EncResult::ERROR;
        }

      /* convert to destination format */
      ret = swr_convert (m_swr_ctx,
                         m_frame->data, dst_nb_samples,
                         (const uint8_t **)frame->data, frame->nb_samples);
      if (ret < 0)
        {
          err = Error ("error while converting");
          return EncResult::ERROR;
        }
      frame = m_frame;

      frame->pts = av_rescale_q (m_samples_count + m_start_pos, (AVRational){1, m_enc->sample_rate}, m_enc->time_base);
      m_samples_count += dst_nb_samples;
    }

  ret = avcodec_send_frame (m_enc, frame);
  if (ret == AVERROR_EOF)
    {
      return EncResult::DONE; // encoder has nothing more to do
    }
  else if (ret < 0)
    {
      err = Error (string_printf ("error encoding audio frame: %s", av_err2str (ret)));
      return EncResult::ERROR;
    }
  for (;;)
    {
      ret = avcodec_receive_packet (m_enc, m_tmp_pkt);
      if (ret == AVERROR (EAGAIN))
        {
          return EncResult::OK; // encoder needs more data to produce something
        }
      else if (ret == AVERROR_EOF)
        {
          return EncResult::DONE;
        }
      else if (ret < 0)
        {
          err = Error (string_printf ("error while encoding audio frame: %s", av_err2str (ret)));
          return EncResult::ERROR;
        }

      /* one packet available */
      if (m_cut_aac_frames)
        {
          m_cut_aac_frames--;
        }
      else if (m_keep_aac_frames)
        {
          ret = write_frame (&m_enc->time_base, m_st, m_tmp_pkt);
          if (ret < 0)
            {
              err = Error (string_printf ("error while writing audio frame: %s", av_err2str (ret)));
              return EncResult::ERROR;
            }
          m_keep_aac_frames--;
        }
    }
}

void
HLSOutputStream::close_stream()
{
  avcodec_free_context (&m_enc);
  av_frame_free (&m_frame);
  av_frame_free (&m_tmp_frame);
  av_packet_free (&m_tmp_pkt);
  swr_free (&m_swr_ctx);
}

Error
HLSOutputStream::open (const string& out_filename, size_t cut_aac_frames, size_t keep_aac_frames, double pts_start, size_t delete_input_start)
{
  assert (m_state == State::NEW);

  avformat_alloc_output_context2 (&m_fmt_ctx, NULL, "mpegts", NULL);
  if (!m_fmt_ctx)
    return Error ("failed to alloc avformat output context");

  /*
   * Since each segment is generated individually, the continuity counter fields of each
   * mpegts segment start at 0, so we expect discontinuities whenever a new segment starts.
   *
   * Players are requested to ignore this by setting this flag.
   */
  int ret = av_opt_set (m_fmt_ctx->priv_data, "mpegts_flags", "+initial_discontinuity", 0);
  if (ret < 0)
    return Error (av_err2str (ret));

  string filename = out_filename;
  if (filename == "-")
    filename = "pipe:1";

  ret = avio_open (&m_fmt_ctx->pb, filename.c_str(), AVIO_FLAG_WRITE);
  if (ret < 0)
    return Error (av_err2str (ret));

  AVDictionary *opt = nullptr;
  const AVCodec *audio_codec;
  Error err = add_stream (&audio_codec, AV_CODEC_ID_AAC);
  if (err)
    return err;

  err = open_audio (audio_codec, opt);
  if (err)
    return err;

  /* Write the stream header, if any. */
  ret = avformat_write_header (m_fmt_ctx, &opt);
  if (ret < 0)
    {
      error ("Error occurred when writing output file: %s\n",  av_err2str(ret));
      return Error ("avformat_write_header failed\n");
    }

  m_delete_input_start = delete_input_start;
  m_cut_aac_frames = cut_aac_frames;
  m_keep_aac_frames = keep_aac_frames;

  // FIXME: correct?
  m_start_pos = pts_start * m_sample_rate - cut_aac_frames * 1024;
  m_start_pos += 1024;

  m_state = State::OPEN;
  return Error::Code::NONE;
}

Error
HLSOutputStream::close()
{
  if (m_state != State::OPEN)
    return Error::Code::NONE;

  // never close twice
  m_state = State::CLOSED;

  Error err;
  while (write_audio_frame (err) == EncResult::OK);
  if (err)
    return err;

  av_write_trailer (m_fmt_ctx);

  close_stream();

  /* Close the output file. */
  if (!(m_fmt_ctx->oformat->flags & AVFMT_NOFILE))
    avio_closep (&m_fmt_ctx->pb);

  /* free the stream */
  avformat_free_context (m_fmt_ctx);

  return Error::Code::NONE;
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

  Error err;
  while (m_audio_buffer.can_read_frames() >= 1024)
    {
      write_audio_frame (err);
      if (err)
        return err;
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
