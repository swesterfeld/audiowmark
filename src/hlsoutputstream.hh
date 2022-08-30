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

#ifndef AUDIOWMARK_HLS_OUTPUT_STREAM_HH
#define AUDIOWMARK_HLS_OUTPUT_STREAM_HH

#include "audiostream.hh"
#include "audiobuffer.hh"

#include <assert.h>

extern "C" {
#include <libavformat/avformat.h>
#include <libavutil/opt.h>
#include <libswresample/swresample.h>
#include <libavutil/avassert.h>
#include <libavutil/timestamp.h>
#include <libavcodec/avcodec.h>
}

class HLSOutputStream : public AudioOutputStream {
  AVStream         *m_st = nullptr;
  AVCodecContext   *m_enc = nullptr;
  AVFormatContext  *m_fmt_ctx = nullptr;

  /* pts of the next frame that will be generated */
  int64_t           m_next_pts = 0;
  int               m_samples_count = 0;
  int               m_start_pos = 0;

  AVFrame          *m_frame = nullptr;
  AVFrame          *m_tmp_frame = nullptr;
  AVPacket         *m_tmp_pkt = nullptr;

  size_t            m_cut_aac_frames = 0;
  size_t            m_keep_aac_frames = 0;

  SwrContext       *m_swr_ctx = nullptr;

  int               m_bit_depth = 0;
  int               m_sample_rate = 0;
  int               m_n_channels = 0;
  AudioBuffer       m_audio_buffer;
  size_t            m_delete_input_start = 0;
  int               m_bit_rate = 0;
  std::string       m_channel_layout;

  enum class State {
    NEW,
    OPEN,
    CLOSED
  };
  State             m_state = State::NEW;

  Error add_stream (const AVCodec **codec, enum AVCodecID codec_id);
  Error open_audio (const AVCodec *codec, AVDictionary *opt_arg);
  AVFrame *get_audio_frame();
  enum class EncResult {
    OK,
    ERROR,
    DONE
  };
  EncResult write_audio_frame (Error& err);
  void close_stream();
  AVFrame *alloc_audio_frame (AVSampleFormat sample_fmt, uint64_t channel_layout, int sample_rate, int nb_samples, Error& err);

  int write_frame (const AVRational *time_base, AVStream *st, AVPacket *pkt);
public:
  HLSOutputStream (int n_channels, int sample_rate, int bit_depth);
  ~HLSOutputStream();

  void set_bit_rate (int bit_rate);
  void set_channel_layout (const std::string& channel_layout);

  Error open (const std::string& output_filename, size_t cut_aac_frames, size_t keep_aac_frames, double pts_start, size_t delete_input_start);
  int bit_depth() const override;
  int sample_rate() const override;
  int n_channels() const override;
  Error write_frames (const std::vector<float>& frames) override;
  Error close() override;
};

#endif /* AUDIOWMARK_HLS_OUTPUT_STREAM_HH */

