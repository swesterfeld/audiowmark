# Audio Watermarking

The **audiowmark** program can add watermarks to audio files and extract
previously embedded watermarks from audio material.
The usage is as follows:


```
usage: audiowmark <command> [ <args>... ]

Commands:
  * create a watermarked wav file with a message
    audiowmark add <input_wav> <watermarked_wav> <message_hex>

  * retrieve message
    audiowmark get <watermarked_wav>

  * compare watermark message with expected message
    audiowmark cmp <watermarked_wav> <message_hex>

  * generate 128-bit watermarking key, to be used with --key option
    audiowmark gen-key <key_file> [ --name <key_name> ]

Global options:
  -q, --quiet             disable information messages
  --strict                treat (minor) problems as errors

Options for get / cmp:
  --detect-speed          detect and correct replay speed difference
  --detect-speed-patient  slower, more accurate speed detection
  --json <file>           write JSON results into file

Options for add / get / cmp:
  --key <file>            load watermarking key from file
  --short <bits>          enable short payload mode
  --strength <s>          set watermark strength              [10]

  --input-format raw      use raw stream as input
  --output-format raw     use raw stream as output
  --format raw            use raw stream as input and output

The options to set the raw stream parameters (such as --raw-rate
or --raw-channels) are documented in the README file.

HLS command help can be displayed using --help-hls
```


# Audiowmark Architecture
<style>	body { max-width: 50em; margin: auto; }	</style>

The **audiowmark** program is used to integrate (`add` command) and extract (`get` command) watermarks (messages of up to 128 bits) into/from audio files.

Internally, the program is organized as nested components, the outermost deals with file IO and command processing.
The commands are implemented via various components that process the watermark, audio signal,
an optional encoding key and user facing information.

## Adding Watermarks

The `audiowmark add <in> <out> <bits> [--key…]` command allows adding watermarks to audio files.
This command takes an audio file, a 128-bit hexadecimal watermark and an optional key as input,
it combines these into a newly generated WAV file. Using the same key, the watermark bits can later
be re-retrieved with the `audiowmark get` command without requiring access to the original
audio input (this is called blind decoding).

By using the encoding key as input, various AES based random number streams are generated to
shuffle, interleave and mix the watermark information into the audio signal.

For robust extraction and forward error correction, the watermark is encoded via convolutional codes
with an order of `15` and a rate of `1/6` (similar to the communication of the
[Mars Pathfinder](https://en.wikipedia.org/wiki/Convolutional_code#Popular_convolutional_codes)).

The expanded watermark bits are transformed into a delta spectrum at a sample rate of 44100Hz and
distributed across various segments (of ca 23 millisecond lengths) of the audio signal and spread
across bands above 800Hz and below 5000Hz.
Based on the delta spectrum, the watermark signal can be modulated and adapted to the current
segment of the input signal before the two are mixed together. To avoid clipping of the output
signal, the final output stage consists of a time local limiter with ca 1 second window.
<!-- TODO: describe the limiter in more detail -->

\pagebreak
An outline of the component interactions to integrate the watermark information via delta
band spectrum into the audio signal is provided in the following chart.

~~~~{.graphviz prog=dot}
digraph "Audiowmark Watermark Embedding" {
  graph[fontsize=13,_fontname="sans"];
  node[fontsize=13,target="_top",_fontname="sans"];
  edge[arrowhead=vee,arrowtail=vee,color="#00000080",_fontname="sans"];
  compound=true;
  concentrate=false;
  rankdir="TB";

bitvec -> get_frame_mod [color=green4];

AudioInput -> Limiter [color=blue];
AudioInput -> in_resampler [color=blue];
AudioInput -> snr_signal_power [style=dashed,color=gold3];

Key -> Random [minlen=2,color=goldenrod];

subgraph cluster_wmadd {
  label=< <b><font face="mono">audiowmark add &lt;AudioInput&gt; &lt;AudioOutput&gt; &lt;Bits&gt; [--key…]                         </font></b> >;

  Random -> get_frame_mod [color=goldenrod];
  in_resampler -> fft_analyzer [color=blue];
  { rank=same Random in_resampler }

  subgraph cluster_WatermarkGen {
    label=< <b>Watermark Generation</b> >; style=dashed; color=grey;
    // fontsize=9; node[fontsize=5,margin=0]; edge[fontsize=5,margin=0];

    fft_analyzer -> apply_frame_mod [color=blue4,xlabel="  513 FFT Bands \r"];
    get_frame_mod -> apply_frame_mod [color=red,xlabel="Up/Down-Band Modulators \r"];
    apply_frame_mod -> wm_synth [color=fuchsia,xlabel="FFT Delta Bands  \r"];

  }

  wm_synth -> out_resampler [color=fuchsia];
  out_resampler -> Limiter [color=fuchsia];
  out_resampler -> snr_signal_power [style=dashed,color=gold3];
}

Limiter -> AudioOutput [color=darkmagenta];
snr_signal_power -> SnrOutput [style=dashed,color=gold3];

bitvec [color=green4,margin=0,label=" Watermark Bits "];
Key [color=goldenrod,margin=0,label=" Key "];
AudioInput [color=blue,margin=0,label=" WAV/MP3 Audio Input File "];
AudioOutput [color=darkmagenta,margin=0,label=" WAV Audio Output File "];
{ rank=same bitvec Key AudioInput }

Random        [color=goldenrod,shape=record,label="Random Stream \l AES128/CTR \l"];
in_resampler  [color=blue,shape=record,label="Resample to 44.1kHz"];
out_resampler [color=fuchsia,shape=record,label="Resample from 44.1kHz"];

Limiter [color=darkmagenta,shape=rect,label=<
<TABLE BORDER="0" ALIGN="LEFT"><TR><TD BALIGN="LEFT">
<b>Mixing ⊕ Limiting</b> <br/>
• Audio and watermark signals are added <br/>
• The result is scaled down to [-0.99…+0.99] <br/>
• Uses 1 second detection window <br/>
</TD></TR></TABLE> >];

snr_signal_power [style=filled,fillcolor="#ffffbb",color=gold3,shape=rect,label=<
<TABLE BORDER="0" CELLBORDER="0" CELLSPACING="1" ALIGN="LEFT"><TR><TD BALIGN="LEFT" CELLPADDING="5">
<b>Power Measurement</b> <br/>
• Signal Power <br/>
• Delta Power <br/>
• Ratio Calculation <br/>
</TD></TR></TABLE> >];
SnrOutput [color=gold3,margin=0,style=filled,fillcolor="#ffffbb",label=" Signal/Noise Ratio Info "];

fft_analyzer [color=blue4,shape=rect,label=<
<TABLE BORDER="0" ALIGN="LEFT"><TR><TD BALIGN="LEFT">
<b>Fourier Transform Analyzer</b> <br/>
• Input: Time domain samples <br/>
• FFT with block Size 1024 <br/>
• Hann Window <br/>
</TD></TR></TABLE> >];

apply_frame_mod [color=fuchsia,shape=rect,label=<
<TABLE BORDER="0" ALIGN="LEFT"><TR><TD BALIGN="LEFT">
<b>Band Modulation ⊗</b> <br/>
• Input bands are Up/Down modulated <br/>
• Factor amounts to ±Amplitude^1% <br/>
• Output: ± Delta bands <br/>
</TD></TR></TABLE> >];

get_frame_mod [color=red,shape=rect,label=<
<TABLE BORDER="0" ALIGN="LEFT"><TR><TD BALIGN="LEFT">
<b>Modulation Frame Generator</b> <br/>
• Encodes Watermark Bits <br/>
• Pregenerate A/B-Blocks <br/>
• Yields A/B-Block frames <br/>
</TD></TR></TABLE> >];

wm_synth [color=fuchsia,shape=rect,label=<
<TABLE BORDER="0" ALIGN="LEFT"><TR><TD BALIGN="LEFT">
<b>Watermark Signal Synthesis</b> <br/>
• Inverse FFT with block size 1024 <br/>
• Cosine window with overlap of 10% <br/>
• Output: Time domain samples <br/>
</TD></TR></TABLE> >];

}
~~~~


At a sample frequency of 44100Hz, the audio signal used for the watermark creation is split into
"Frames" of 1024 samples each, which corresponds to segments of ca 23 millisecond length. These
frames are transformed from the time domain (samples) into the frequency domain (spectral bands)
and vice versa to apply the watermark embedding in certain spectral bands.

Data and synchronization bits are encoded across several frames, with different levels of
redundancy. In the "Modulation Frame Generator" the number of frames that compose all encoded
information needed to find and extract the watermark bits are combined into two types of "Blocks".

\pagebreak
A detailed chart of the component interactions for the Frame and Block generation in the
"Modulation Frame Generator" is provided in the next chart.

~~~~{.graphviz prog=dot}
digraph "Modulation Frame Generator" {
  graph[fontsize=13,_fontname="sans"];
  node[fontsize=13,target="_top",_fontname="sans"];
  edge[arrowhead=vee,arrowtail=vee,color="#00000080",_fontname="sans"];
  compound=true;
  concentrate=false;
  rankdir="TB";

Random -> ab_generators [color=goldenrod,lhead=cluster_ModulationFrameGenerator];
//Random -> randomize_bit_order [color=goldenrod,xlabel="bit_order R5"];
//Random -> mark_sync [color=goldenrod,xlabel="sync_up_down R2"];
//Random -> mark_data [color=goldenrod,xlabel="data_up_down R1"];
//Random -> frame_pos [color=goldenrod,xlabel="frame_position R6"];

bitvec -> conv_encode [color=green4,minlen=2];
{ rank=same Random bitvec }

subgraph cluster_ModulationFrameGenerator {
  label=< <b>Modulation Frame Generator</b> >; style=dashed; color=red;

  ab_generators -> conv_encode [color=green4,minlen=2];
  conv_encode -> randomize_bit_order [color=green4,minlen=2];
  { rank=same ab_generators conv_encode }

  UpDownGen -> mark_sync [color=goldenrod];
  UpDownGen -> gen_mix_entries [color=goldenrod];
  frame_pos -> mark_sync [color=goldenrod,minlen=1];
  // --linear: frame_pos -> mark_data;
  frame_pos -> gen_mix_entries [color=goldenrod];

  randomize_bit_order -> mark_data [color=goldenrod];
  init_frame_mod_vec -> get_frame_mod [color=cyan4];
  mark_sync -> init_frame_mod_vec [color=cyan3,minlen=1];
  mark_data -> init_frame_mod_vec [color=teal];
  gen_mix_entries -> mark_data [color=goldenrod];
  { rank=same mark_data mark_sync }
}

get_frame_mod -> apply_frame_mod [color=red,xlabel="Frame Up/Down-Band Modulators",minlen=2];
apply_frame_mod [shape=plain,label=" "];

bitvec [color=green4,margin=0,label=" Watermark Bits "];
Random [color=goldenrod,shape=record,label="
Random Stream \l
AES128/CTR \l
Streams R1…R6 \l
"];

conv_encode [color=green4,shape=rect,label=<
<TABLE BORDER="0" ALIGN="LEFT"><TR><TD BALIGN="LEFT">
<b>Convolutional Code Expansion</b> <br/>
• Pads watermark with termination zeros <br/>
• Combines bit stream with A/B constants <br/>
• Generates output stream of parity bits <br/>
• Generates 858 parity bits A-Block <br/>
• Generates 858 parity bits B-Block <br/>
</TD></TR></TABLE>>];

ab_generators [color=green4,shape=rect,label=<
<TABLE BORDER="0" ALIGN="LEFT"><TR><TD BALIGN="LEFT">
<b>Convolutional Code Parameters</b> <br/>
• Convolutional code with rate 1/6 <br/>
• Order 15, needs 15 termination bits <br/>
• Six constants for A-Block and B-Block <br/>
• Forward correction of ca ≈20% bit errors <br/>
• Encodes 128 bits in 858 bit blocks <br/>
</TD></TR></TABLE>>];

UpDownGen [color=goldenrod,shape=rect,label=<
<TABLE BORDER="0" ALIGN="LEFT"><TR><TD BALIGN="LEFT">
<b>Up/Down-Band Generator</b> <br/>
• Uses per-frame shuffling seed <br/>
• Picks random bands, 30 UP, 30 DOWN <br/>
• Bands are between ca 861Hz…4307Hz <br/>
</TD></TR></TABLE>>];

mark_sync [color=cyan3,shape=rect,label=<
<TABLE BORDER="0" ALIGN="LEFT"><TR><TD BALIGN="LEFT">
<b>Synchronization Frame Generator</b> <br/>
• Encodes 6 sync bits in 6 * 85 frames <br/>
• A-Block bit pattern: 010101 <br/>
• B-Block bit pattern: 101010 <br/>
• Randomizes Up/Down-Band shifts [R2] <br/>
• Output: 510 Frames * 60 Up/Down-Bands <br/>
</TD></TR></TABLE>>];

gen_mix_entries [color=goldenrod,shape=rect,label=<
<TABLE BORDER="0" ALIGN="LEFT"><TR><TD BALIGN="LEFT">
<b>Mix Entry Generator (skipped for --linear)</b> <br/>
• Generates list of data bit encoding bands <br/>
• Uses 30 up + 30 down bands in 2 frames per bit <br/>
• Randomizes Up/Down-Band shifts [R1] <br/>
• Shuffles data bit association of entries [R4] <br/>
• Output: 2 * 858 * 30 Up/Down band pairs <br/>
</TD></TR></TABLE>>];

mark_data [color=teal,shape=rect,label=<
<TABLE BORDER="0" ALIGN="LEFT"><TR><TD BALIGN="LEFT">
<b>Data Frame Generator</b> <br/>
• Encodes 858 data bits in 858 * 2 frames <br/>
• Encodes A-Blocks, B-Blocks in turn <br/>
• Omits Mix Entry Generator with --linear <br/>
• Randomizes Up/Down-Band shifts [R1] <br/>
• Output: 1716 Frames * 60 Up/Down-Bands <br/>
</TD></TR></TABLE>>];

frame_pos [color=goldenrod,shape=rect,label=<
<TABLE BORDER="0" ALIGN="LEFT"><TR><TD BALIGN="LEFT">
<b>Frame Position Randomization</b> <br/>
• Mixes sync + data frames <br/>
• Shuffles frame positions <br/>
• Uses random stream [R6] <br/>
</TD></TR></TABLE>>];

randomize_bit_order [color=goldenrod,shape=rect,label=<
<TABLE BORDER="0" ALIGN="LEFT"><TR><TD BALIGN="LEFT">
<b>Randomize Bit Order for ENCODE</b> <br/>
• Reversible shuffle for encode/decode <br/>
• Shuffles/interleaves bit stream [R5] <br/>
• Interleaving improves robustness <br/>
• Reduces bit stream damage impact <br/>
</TD></TR></TABLE>>];

init_frame_mod_vec [color=cyan4,shape=rect,label=<
<TABLE BORDER="0" ALIGN="LEFT"><TR><TD BALIGN="LEFT">
<b>A/B-Block Frame Modulator Composition</b> <br/>
• Interleaves synchronization and data frames <br/>
• Pulls and interleaves each block type separately <br/>
• Output: Up/down band modulators for 1 block <br/>
</TD></TR></TABLE>>];

get_frame_mod [color=red,shape=rect,label=<
<TABLE BORDER="0" ALIGN="LEFT"><TR><TD BALIGN="LEFT">
<b>Modulation Frame Selector</b> <br/>
• Yields A-Block band modulators per frame <br/>
• Yields B-Block band modulators and starts over <br/>
• Output: Up/down band modulators for 1 frame <br/>
</TD></TR></TABLE>>];

}
~~~~

The watermark is encoded and embedded into the audio signal in two block types, A-Blocks and B-Blocks.
Given ideal transmissions, the watermark can be extracted from each of the block types.
In case of distorted and noisy transmissions where watermark extraction from either block type fails,
a combination of segments with A-Block and B-Block data may still lead to successful recovery of the original watermark.

In order to support watermark extraction from clipped excerpts of the input stream, a fixed pattern
of synchronization bits is integrated into the data blocks with much higher redundancy than the data bits.
The fixed pattern allows detection of A-Blocks and B-Blocks as such to aid the watermark extraction.

The user provided encoding `Key` seeds an AES based pseudo random number generator in
[Counter Mode](https://en.wikipedia.org/wiki/Block_cipher_mode_of_operation#CTR)
that is used to determine encoding places, randomize the noise introduced by the watermark and to interleave encoding
for robustness. Without the key, the watermark information cannot be retrieved and its presence can not be detected.
The different types of random data streams used for the distribution of the embedded watermark information are as follows:

* R1 - Used to randomizes Up/Down band shifts for watermark data bits.
* R2 - Used to randomizes Up/Down band shifts for watermark synchronization bits.
* R3 - Currently unused.
* R4 - Used to mix (shuffle) data bit associations of Up/Down bands distributed across several frames.
* R5 - Used to shuffle (interleave) the bit stream. Due to redundancy in the generated bit stream, interleaving reduces the
  number of adversely affected bits by bursts (holes) in transmission loss.
* R6 - Used to randomize and mix data frames with synchronization frames, this makes synchronization frames unlikely to be detectable without the encoded key.

## Extracting Watermarks

The `audiowmark get <watermarked_wav> [--key…]` command extracts a watermark from an audio file.
This command takes an audio file and an optional key as input.
With the same key used during watermark embedding, synchronization bits are determined and searched
for in the audio content.
If synchronization bit matches are detected, encoded watermark information can be located,
extracted and decrypted with error correction.
The retrieval does not require access to the original audio input (this is called blind decoding).
The detection results are produced on *stdout* with accompanying information about the location,
match quality and a measure for likely decoding errors.

An outline of the component interactions to locate and extract the watermark information from the
frequency spectrum in the audio signal is provided in the following charts.

~~~~{.graphviz prog=dot}
digraph "Audiowmark Watermark Extraction" {
  newrank=true;
  graph[fontsize=13,_fontname="sans"];
  node[fontsize=13,target="_top",_fontname="sans"];
  edge[arrowhead=vee,arrowtail=vee,color="#00000080",_fontname="sans"];
  compound=true;
  concentrate=false;
  rankdir="TB";

WavData -> in_resampler [color=darkmagenta];

Key -> Random [minlen=2,color=goldenrod];

subgraph cluster_wmget {
  label=< <b><font face="mono">audiowmark get &lt;AudioFile&gt; [--key…]
                                                    </font></b> >;

  { rank=same Random in_resampler }
  in_resampler -> fft_range [color=darkmagenta];
  Random -> conv_decode_soft [minlen=2,color=goldenrod,
    lhead=cluster_WatermarkExtraction]; // fake arrow target for cluster alignment

  subgraph cluster_WatermarkExtraction {
    label=< <b>Watermark Extraction</b> >; style=dashed; color=green4;

    BlockDecoder -> ClipDecoder [color=darkorchid3];
    fft_range -> BlockDecoder [color=blue4,constraint=1];

    { rank=same fft_range conv_decode_soft }
    conv_decode_soft -> BlockDecoder [color=green4,dir=both];
    conv_decode_soft -> ClipDecoder [color=green4,dir=both];
  }
  BlockDecoder -> SyncFinder [color=cyan4,minlen=4,dir=both,xlabel="Mode::BLOCK"];
  ClipDecoder -> SyncFinder  [color=cyan4,minlen=1,dir=both,xlabel="Mode::CLIP",constraint=false];

  subgraph cluster_SyncFinder {
    label=""; style=dashed; color=cyan4; node[margin=0]; edge[margin=0]; margin=0;
    SyncFinder [shape=plaintext,
    label=< <table border="0"><tr><td BALIGN="LEFT">
<b>Synchronization Position Finder </b> <br/>
• Performs coarse search for synchronization <br/>
   bit markers in the frequency spectrum <br/>
• Searches for A/B-Block synchronizations <br/>
• Calculates score for possible locations <br/>
   and picks the 5 best matches <br/>
• Refines the exact block locations with a fine <br/>
   grained search for synchronization markers <br/>
• Short audio segments are dealt with by <br/>
   adding zero padding (in Mode::CLIP) <br/>
</td></tr></table> >];
  }


  ClipDecoder -> result_set_print [color=darkorchid4];

  { rank=same SyncFinder result_set_print }
}


result_set_print -> bitvec [color=green4];

Key [color=goldenrod,margin=0,label=" Key "];
WavData [color=darkmagenta,margin=0,label=" WAV/MP3 Audio Input File "];
{ rank=same Key WavData }
in_resampler  [color=darkmagenta,shape=record,label="Resample to 44.1kHz"];
Random        [color=goldenrod,shape=record,label="Random Stream \l AES128/CTR \l"];

fft_range [color=blue4,shape=rect,label=< <table border="0" align="left"><tr><td balign="left">
<b>Fourier Transform Analyzer</b> <br/>
• Input: Frames with time domain samples <br/>
• Short inputs are zero-padded as needed (Mode::CLIP) <br/>
• FFT with block Size 1024 <br/>
• Hann Window <br/>
• Output: 513 FFT Bands <br/>
</td></tr></table> >];

BlockDecoder [color=darkorchid3,shape=rect,label=< <table border="0" align="left"><tr><td balign="left">
<b>Best Block Decoder</b> <br/>
• Detect A/B-Blocks via 'Synchronization Position Finder' in Mode::BLOCK <br/>
• Only works for audio clips with at least 52 seconds and proper block <br/>
   alignment (or up from 104 seconds without alignment)  <!-- 1024 * (510 + 1716) / 44100 --> <br/>
• Reconstruct Up/Down-Band associations (reverses 'Mix Entry Generator') [R4] <br/>
• Estimate bit vectors resulting from average deviations in <br/>
   randomized Up/Down-Band shifts [R1] <!-- mix_decode --> <br/>
• Unshuffle bit order (reverses 'Randomize Bit Order for ENCODE') [R5]<br/>
• Normalize soft bits (normalization of estimated bit vectors) <br/>
• Reconstruct watermark bits using 'Soft-Decision Decoder' <br/>
• Decode individual A-Blocks, B-Blocks and if present AB-Blocks <br/>
• As last resort, attempt an AB-Block decode on accumulated bit vectors <br/>
   averaged over all selected blocks
</td></tr></table> >];

conv_decode_soft [color=green4,shape=rect,label=< <table border="0" align="left"><tr><td balign="left">
<b>Soft-Decision Decoder</b> <br/>
• Utilize Viterbi algorithm <br/>
• Decode blocks of 858 parity bits <br/>
• Reconstructs 128 payload bits <br/>
• Uses 'Convolutional Code Parameters' <br/>
• Decode A-, B- and AB-Blocks <br/>
</td></tr></table> >];

ClipDecoder [color=darkorchid4,shape=rect,label=< <table border="0" align="left"><tr><td balign="left">
<b>Short Audio Clip Decoder</b> <br/>
• Used only for small audio clips up to 160 seconds (3.1 blocks) <!-- 1024 * (510 + 1716) / 44100 * 3.1 --> <br/>
• Uses zero padding (silence) around short audio clips to construct 3 blocks <br/>
• Zero padded regions are ignored during synchronization detection and scoring <br/>
• Determine alignment with 'Synchronization Position Finder' in Mode::CLIP <br/>
• Select the five blocks with the best detected synchronization markers <br/>
• Reconstruct Up/Down-Band associations (reverses 'Mix Entry Generator') [R4] <br/>
• Estimate bit vectors resulting from average deviations in <br/>
   randomized Up/Down-Band shifts [R1] <!-- mix_decode --> <br/>
• Unshuffle bit order (reverses 'Randomize Bit Order for ENCODE') [R5]<br/>
• Normalize soft bits (normalization of estimated bit vectors) <br/>
• Reconstruct watermark bits using 'Soft-Decision Decoder' <br/>
• Attempt AB-Block decoding at alignment position <br/>
• Attempt secondary AB-Block decoding in case of audio data surplus <br/>
• In practice, ca 10 seconds are needed for reliable detection, <br/>
   in good scenarios up to 3 seconds may suffice <br/>
</td></tr></table> >];


result_set_print [style=filled,fillcolor="#eeffdd",color=green4,shape=rect,label=<
<TABLE BORDER="0" ALIGN="LEFT"><TR><TD BALIGN="LEFT">
<b>Result Set Printing</b> <br/>
• Sort results by block type and score <br/>
• Print potential watermarks on stdout <br/>
• Provide time offset, score quality, <br/>
   block type and a measure for likely <br/>
   decoding errors per watermark <br/>
</TD></TR></TABLE> >];

bitvec [style=filled,fillcolor="#eeffdd",color=green4,margin=0,label=" Watermark Bits "];

}
~~~~

At a sample frequency of 44100Hz, a spectral analysis is performed on the audio signal and the spectrum is then searched for known synchronization markers.
Upon detection of A/B-Block synchronization positions, watermark bits are extracted from known data bit locations, while making use of the embedded redundancy to make the detection more robust.
Due to high redundancy and wide spread of watermark information, bits often can still be extracted from audio clips that are heavily shortened. To employ the full detection machinery to very short clips, symmetric zero padding is used to provide enough input samples (zero padded regions are ignored during scoring however).

Since detection success is directly dependent on the precise bit stream synchronization, an iterative process is used for fast approximation of synchronization locations with later refinements to yield precise results.

The purpose of the synchronization algorithm is to find the location of the
watermark A/B blocks in the input signal. This is important because the
signal may have been cropped so that the location of the blocks is not
known. To be able to find the locations of the blocks, while adding the
watermark, some sync bits are added to each block with relatively high redundancy.
The values of these sync bits are known, for an A block they are 010101, for a
B block they are 101010. The up- and down-bands used for the sync bits and
offsets of all frames that belong to sync bits inside the A / B block are known
and determined by the key.

To perform the actual synchronization and locate the start of an A (or B) block,
two steps are performed.

* As a first step, the synchronization algorithm tests all possible positions
for the start of an A (or B) block using a step size of 256 (1/4 frame size)
and tries to decode the sync bits at the expected locations relative to the
start of the block. Since the values and locations of the bits are known, a sync
score can be computed that indicates how good the bits in the actual audio
input at this position match the expected bit sequence.

* For all start locations with a significantly high sync score, in a second
step the actual start position is searched by trying all different start
locations near to the original match with a smaller finer step size. Again
a sync score can be computed and compared to a second threshold to decide
whether this is location is really likely to contain a data block. If the match
is good enough the start location will be used to decode the data bits in the
block.

Besides using this strategy to find "whole" data blocks, there is also a
variant of the synchronization algorithm that is used if the audio signal
is very short. It can find the location of the watermark even if the length
of the input signal is too short to contain a complete data block. To be
able to do this, the input signal is zero padded before sync detection and
then the usual algorithm to find whole blocks is used.

The following chart provides the detail of the steps involved in determining the synchronization locations.


~~~~{.graphviz prog=dot}
digraph "Synchronization Position Finder" {
  graph[fontsize=13,_fontname="sans"];
  node[fontsize=13,target="_top",_fontname="sans"];
  edge[arrowhead=vee,arrowtail=vee,color="#00000080",_fontname="sans"];
  compound=true;
  concentrate=false;
  rankdir="TB";

  in_resampler;
  Random;
  { rank=same in_resampler Random }

  Random -> frame_pos [color=goldenrod,lhead=cluster_SyncFinder,minlen=2];

  subgraph cluster_SyncFinder {
    label=< <b>Synchronization Position Finder</b> >; style=dashed; color=cyan4;

    in_resampler -> fft_analyzer [color=darkmagenta];
    fft_analyzer -> sync_fft_256_8 [color=blue4,dir=both];
    frame_pos -> init_up_down [color=goldenrod];
    UpDownGen -> init_up_down [color=goldenrod];
    init_up_down -> search_approx [color=cyan3];
    search_approx -> sync_select_by_threshold [color=cyan3];
    sync_select_by_threshold -> search_refine [color=cyan3];

    sync_fft_256_8 -> search_refine [style=dashed,dir=back,color=darkcyan,xlabel=" Refining\l Feedback\l                "];
    sync_fft_256_8 -> sync_decode [color=darkcyan,style=dashed,label="Repeat\lRefined\l                "];
    sync_decode -> search_refine [color=darkcyan,style=dashed,label="Repeat\lRefined\l                "];

    sync_fft_256_8 -> sync_decode [color=cyan3];
    sync_decode -> search_approx [color=cyan3];

    { rank=same fft_analyzer init_up_down }
  }

search_refine -> sync_scores [color=cyan4,xlabel="Scoring and A/B-Type for potential blocks",minlen=2];
sync_scores [shape=plain,label=" "];

in_resampler  [color=darkmagenta,shape=record,label="
WAV/MP3 Audio Input File \l
\l
Resampled to 44.1kHz \l
"];

Random [color=goldenrod,shape=record,label="
Random Stream \l
AES128/CTR \l
Streams R1…R6 \l
"];
{ rank=same in_resampler Random }

frame_pos [color=goldenrod,shape=rect,label=<
<TABLE BORDER="0" ALIGN="LEFT"><TR><TD BALIGN="LEFT">
<b>Frame Position Randomization</b> <br/>
• Mixes sync + data frames <br/>
• Shuffles frame positions <br/>
• Uses random stream [R6] <br/>
</TD></TR></TABLE> >];

UpDownGen [color=goldenrod,shape=rect,label=<
<TABLE BORDER="0" ALIGN="LEFT"><TR><TD BALIGN="LEFT">
<b>Up/Down-Band Generator</b> <br/>
• Uses per-frame shuffling seed <br/>
• Picks random bands, 30 UP, 30 DOWN <br/>
• Bands are between ca 861Hz…4307Hz <br/>
</TD></TR></TABLE> >];

init_up_down [color=cyan3,shape=rect,label=<
<TABLE BORDER="0" ALIGN="LEFT"><TR><TD BALIGN="LEFT">
<b>Synchronization Bit Frame Generator</b> <br />
<b>(Mode::BLOCK &amp; Mode::CLIP)</b> <br/>
• Generates 6 sync bits in 6 * 85 frames (* 2 for Mode::CLIP) <br/>
• Randomizes Up/Down-Band shifts [R2] <br/>
• Sorts synchronization bit frames by frame index <br/>
• Mode::BLOCK Output: 510 Bit Frames with 60 Up &amp; 60 Down-Bands <br/>
• Mode::CLIP Output: 1020 Bit Frames with 60 Up &amp; 60 Down-Bands <br/>
</TD></TR></TABLE> >];

fft_analyzer [color=blue4,shape=rect,label=<
<TABLE BORDER="0" ALIGN="LEFT"><TR><TD BALIGN="LEFT">
<b>Fourier Transform Analyzer</b> <br/>
• Input: Time domain samples <br/>
• FFT with block Size 1024 <br/>
• Hann Window <br/>
• Output: 513 FFT Bands <br/>
</TD></TR></TABLE> >];

sync_fft_256_8 [color=blue4,shape=rect,label=<
<TABLE BORDER="0" ALIGN="LEFT"><TR><TD BALIGN="LEFT">
<b>Decibel Quantifier</b> <br/>
• Uses coarse stepping of 256 values for approximate search <br/>
• A stepping of 256 equates 1/4th FFT block <br/>
• Uses fine stepping of 8 values for refined search <br/>
• Pulls FFT Bands for all input blocks <br/>
• Computes dB for all bands of all blocks <br/>
</TD></TR></TABLE> >];

sync_decode [color=cyan3,shape=rect,label=<
<TABLE BORDER="0" ALIGN="LEFT"><TR><TD BALIGN="LEFT">
<b>Synchronization Bit Matching</b> <br/>
• Collect Up/Down-Band magnitudes for sync bits <br/>
• Determine match quality for alternating bit patterns <br/>
• Apply watermark strength dependent thresholds <br/>
• Decide A/B-Block based on 010101 / 101010 detection <br/>
</TD></TR></TABLE> >];

search_approx [color=cyan3,shape=rect,label=<
<TABLE BORDER="0" ALIGN="LEFT"><TR><TD BALIGN="LEFT">
<b>Approximate Synchronization Frame Search</b> <br/>
• Skips over zero-padded samples (Mode::CLIP) <br/>
• Computes multiple time-shifted FFT vectors <br/>
• Uses coarse subframe stepping of 256 values <br/>
• Overlaps frames for sync detection by 1/4th frame <br/>
• Scores positions for synchronization matches <br/>
</TD></TR></TABLE> >];

sync_select_by_threshold [color=cyan3,shape=rect,label=<
<TABLE BORDER="0" ALIGN="LEFT"><TR><TD BALIGN="LEFT">
<b>Synchronization Frame Selection</b> <br/>
• Due to the subframe stepping, good and bad matches can be expected to alternate <br/>
• Identification of local match maxima <br/>
• Strength dependent threshold determines minimum match quality <br/>
• Selection of likely match positions via maxima and threshold <br/>
</TD></TR></TABLE> >];

search_refine [color=darkcyan,shape=rect,label=<
<TABLE BORDER="0" ALIGN="LEFT"><TR><TD BALIGN="LEFT">
<b>Refined Synchronization Frame Search</b> <br/>
• Computes fine-stepped time-shifted FFT vectors around selected frames <br/>
• Searches ±16 subframes around previously detected good scores <br/>
• Keeps subframe if the score (synchronization frame detection quality) improves <br/>
</TD></TR></TABLE> >];

}
~~~~

The user provided 128-Bit AES key is essential to determine spectral bands, encoding patterns, and bit locations.
During decoding, the same Pseudo Random Number Generator sequences R1…R6 are used that facilitated watermark embedding.
By using the same AES key and a cryptographically secure PRNG, the sequences are uniformly distributed and deterministically reproducible but cannot be extrapolated.
This prevents watermark extraction or modification by anyone without possession of the exact encoding key.

## The Patchwork Algorithm

![Example Spectrum](example-spectrum.png)

To store one single bit inside a spectrum, **audiowmark** uses the patchwork
algorithm. From the frequency bands of the spectrum (generated by computing the
FFT of one frame), two groups are choosen in the frequency range of the watermark
using the pseudo random number generator. These are called up- and down-bands.
In the example above, the up-bands are red and the down-bands are green.
Typically there are 30 up- and 30 down-bands and the other bands do not carry
information.

To embed a single bit, the following changes are made to the spectrum:

 * to **store a 1 bit**, each magnitude of each up-band is increased by a small amount,
   and each magnitude of each down-band is decreased by a small amount (this is
   shown by the small arrows in the example image)

 * to **store a 0 bit**, each magnitude of each up-band is decreased, and each magnitude of
   each up-band is increased (the opposite of the small arrows in the example image)

Since we have pseudo-randomly choosen the up- and down-bands from the spectrum,
we can expect that if we sum up all values of the up-bands and sum up all
values of the down-bands **before** embedding the bit, we will get a similar
result (because the mean value of all spectrum bins is shared between the two).

However, since we increased all elements of the up-bands and decreased all
elements of the down-bands **after embedding a 1 bit**, the sum of the up-bands
should be **greater than** the sum of the down-bands.

So to decode the bit from the spectrum, we can simply use the rule

 * **decode as 1 bit**, if the sum of the up-bands is greater than the sum
   of the down-bands

 * **decode as 0 bit**, if the sum of the up-bands is smaller than the sum
   of the down-bands

In the actual implementation, increasing/decreasing the magnitude of the
up-/down-bands is done by generating a watermark signal with the right
magnitude/phase for each frame that only contains the changes. So we
compute a delta spectrum, which is then passed to the IFFT, windowed and then
added to the original audio, so that the sum has the desired modified spectrum
magnitude.

The detection is performed on dB values of the magnitudes of the spectrum
obtained from the FFT, so the sums of the dB values of up-/down-bands are
computed and compared to decide whether a 0 bit or 1 bit was received.

The patchwork algorithm does not guarantee that encoding/decoding will always
yield the right result at the lowest level of embedding/decoding one bit (as
the difference of the up-/down-bands can be too big before embedding due to
the original signal). However error correction and redundancy by embedding a
bit in more than one frame makes the whole process reliable at a higher level.

There are three improvements over the basic patch work algorithm described
above, which make the watermark detection more accurate:

* To use soft-decoding for the convolutional decoder, instead of deciding
whether a 0 or 1 bit was received by comparing the two sums directly before
decoding the convolutional code to obtain the message bits, the difference
between the two sums is normalized and is used as a soft-bit input for the
Viterbi algorithm.

* Instead of storing one data bit in each frame spectrum, a data bit uses up-
and down-bands from different frames. This is called mix-encoding, which
spreads the information of each data bit over many frames.

* As described above, the original signal can have some negative effect
on the performance of the decoder, since the sum of the up-bands and the
sum of the down-bands will be different even before embedding the bits.
To make detection more reliable, the original signal level for each bin is
estimated by taking the average value of the previous and next spectrum and
subtracted before computing the sum of the up- and down-bands.

## Speed Detection

As one of the later developments, a dedicated speed detection facility has been integrated that explores the ability to extract watermarks from
audio segments played back at unknown rates.
In scenarios where audio has been resampled and pitched at a constant rate, synchronization markers may still be
detectable by searching the audio content at varying resampled playback rates.

The `audiowmark --detect-speed` command line option attempts to detect playback rate changes compared to the original material used for embedding within 80% to 125%.

Detection of playback rate modifications is approached in several steps.
First, the detector picks two short audio clips (ca 25 seconds) with high signal energy and performs multiple coarse scans while detecting <0.2% rate modulations.
This rough speed estimate is improved upon with secondary scans around 1/20‰ rate modulations on an audio clip of 50 second.

Executing coarse and fine detection runs at varying resampled playback rates with multiple refining steps consumes a lot of processing resources.
To speed up the detection, resampling and scanning runs are carried out on a downsampled version of the audio material (by a factor of 2) and detection runs are parallelized across all available CPU cores.

Finally, the watermark extraction is carried out on a resampled version of the audio material at the most likely detected playback rate in addition to regular watermark detection, because the detected playback rate may have been guessed wrongly.


<!--
# https://graphviz.org/documentation/
BUILD:
apt install -y python3-pygraphviz python3-pandocfilters
pandoc -F graphviz.py audiowmark.md -o audiowmark.html
pandoc -F graphviz.py audiowmark.md -V papersize:a4 -V geometry:margin=2cm -o audiowmark.pdf
-->
