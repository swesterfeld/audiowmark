% VIDEOWMARK(1) videowmark | User Commands

# NAME

videowmark — video watermarking tool (using audio track)

# SYNOPSIS

**videowmark** [*options*] *command* [*args...*]

# DESCRIPTION

**videowmark** creates watermarked video files containing an embedded message,
and can retrieve an embedded message from a previously watermarked video.

Messages are provided as hexadecimal data and embedded in the audio track of
the video, intended to be robust against common audio processing operations.

# COMMANDS

**add** *input_video* *watermarked_video* *message_hex*

: Create a watermarked video file containing the provided message.

**get** *watermarked_video*

: Retrieve the embedded message from a watermarked video and print it to
standard output.

# GLOBAL OPTIONS

These options apply to all commands.

**\--strength** *s*

: Set watermark strength.

**\--key** *file*

: Load watermarking key material from *file*.

**-q**, **\--quiet**

: Disable informational messages.

**-v**, **\--verbose**

: Enable verbose output from the underlying **ffmpeg** processing.

# ARGUMENTS

*input_video*

:   Path to the input (unwatermarked) video file.

*watermarked_video*

:   Path to the output (watermarked) video file (for **add**), or the input
    watermarked video file (for **get**).

*message_hex*

:   Message to embed, provided as hexadecimal bytes.

# EXAMPLES

Embed a message into a video:

`videowmark add input.mp4 output.mp4 68656c6c6f`

Retrieve a message from a watermarked video:

`videowmark get output.mp4`

# EXIT STATUS

**0**

:   The operation completed successfully.

**>0**

:   An error occurred.

# SEE ALSO

- **ffmpeg**(1)
- **audiowmark**(1)

# HOMEPAGE

https://uplex.de/audiowmark/

# AUTHOR

**videowmark** is part of the audiowmark project and was written by
Stefan Westerfeld <stefan@space.twc.de>.

# COPYRIGHT

Copyright © 2018–2020 Stefan Westerfeld <stefan@space.twc.de>.
Licensed under the GNU General Public License, version 3 or later.
