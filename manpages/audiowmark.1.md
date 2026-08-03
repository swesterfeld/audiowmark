% AUDIOWMARK(1) audiowmark | User Commands

# NAME

audiowmark — audio watermarking tool

# SYNOPSIS

**audiowmark** [*options*] *command* [*args...*]

# DESCRIPTION

**audiowmark** is an open-source C++ audio watermarking utility that embeds and retrieves
128-bit messages in audio files in a way that is inaudible to most listeners.

It supports blind decoding (watermark extraction without access to the
original file) and robust forward error correction, allowing watermarks to
survive lossy compression formats such as MP3 and OGG.

# COMMANDS

**add** *input_wav* *watermarked_wav* *message_hex*

: Create a watermarked WAV file containing the provided message.

**get** *watermarked_wav*

: Retrieve the embedded watermark message from a watermarked WAV file and print
it to standard output.

**cmp** *watermarked_wav* *message_hex*

: Compare the embedded watermark message against the expected message.

**gen-key** *key_file* [ **\--name** *key_name* ]

: Generate a 128-bit watermarking key and write it to *key_file*.

# GLOBAL OPTIONS

These options apply to all commands.

**-q**, **\--quiet**

: Disable informational messages.

**\--strict**

: Treat minor problems as errors.

# OPTIONS FOR GET / CMP

**\--detect-speed**

: Detect and correct replay speed differences.

**\--detect-speed-patient**

: Slower but more accurate replay speed detection.

**\--json** *file*

: Write JSON-formatted results to *file*.

# OPTIONS FOR ADD / GET / CMP

**\--key** *file*

: Load watermarking key material from *file*.

**\--short** *bits*

: Enable short payload mode using the specified number of bits.

**\--strength** *s*

: Set watermark strength.
The default value is **10**.

**\--input-format raw**

: Use a raw audio stream as input.

**\--output-format raw**

: Use a raw audio stream as output.

**\--format raw**

: Use raw audio streams for both input and output.

Options to configure raw stream parameters (such as sample rate or channel
count) are documented in the README.

# EXAMPLES

Embed a watermark into an audio file:

`audiowmark add input.wav output.wav 00112233445566778899aabbccddeeff`

Extract a watermark from an audio file:

`audiowmark get output.wav`

Compare a watermark with an expected message:

`audiowmark cmp output.wav 00112233445566778899aabbccddeeff`

Generate a watermarking key:

`audiowmark gen-key example.key --name "example key"`

# EXIT STATUS

**0**

:   The operation completed successfully.

**>0**

:   An error occurred.

# SEE ALSO

- **ffmpeg**(1)
- **sox**(1)

# HOMEPAGE

https://uplex.de/audiowmark/

# AUTHOR

audiowmark was written by Stefan Westerfeld <stefan@space.twc.de>.

# COPYRIGHT

Copyright © 2018–2020 Stefan Westerfeld <stefan@space.twc.de>.
Licensed under the GNU General Public License, version 3 or later.
