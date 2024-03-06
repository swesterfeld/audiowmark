/*
 * Copyright (C) Andreas Strohmeier
 *
 * This program is a port of the Linux bash script videowmark,
 * which is part of the audiowmark program by Stefan Westerfeld,
 * into the C++ programming language.
 * To keep it simple, there is only this single CPP file without a
 * header file.
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
 *
 */

 //---------------------------------------------------------------------------

#include <iostream>
#include <string>
#include <cstdio>
#include <windows.h>
#include <algorithm>
#include <filesystem>

//---------------------------------------------------------------------------

using namespace std;

//---------------------------------------------------------------------------

string g_sVersion = "videowmark 0.0.5";
string g_sFFMPEG_VERBOSE = "-v error";
string g_sFFProbe = "ffprobe.exe";
string g_sFFMpeg = "ffmpeg.exe";
string g_sAudiowmark = "audiowmark.exe";
int g_iQuiet = 0;

//---------------------------------------------------------------------------

#define STRING_LENGTH 4096

//---------------------------------------------------------------------------

//show message and exit
int die(string sErrorMsg)
{
    printf("videowmark: error: %s\n", sErrorMsg.c_str());
    exit(1);
}
//---------------------------------------------------------------------------

// Fix me
// Maybe there is es better solution
//
// The result of getenv is not correct in cygwin environment.
//
// Example: getenv("TEMP");
//
// Result : /cygdrive/c/path/to/temp
// But should be: c:\path\to\temp
//
// This is a workaround to repair the path
string repair_cygwin_path(string sPath)
{
int i = 0;
string sResult = "";

    if(int(sPath.find("/cygdrive/")) == 0 )
    {
        sPath.erase(0, 10);

        for(i = 0; i < int(sPath.length()); i++)
        {
            if(sPath[i] == '/')
            {
                sPath[i] = '\\';
            }
        }

        if(sPath[1] != ':')
        {
            sPath.insert(1, ":");
        }
    }

    sResult = sPath;

    return sResult;
}
//---------------------------------------------------------------------------

// Converts the Windows path to the UNIX path
string get_unix_path(string sPath)
{
int i = 0;
string sResult = "";

    for(i = 0; i < int(sPath.length()); i++)
    {
        if(sPath[i] == '\\')
        {
            sPath[i] = '/';
        }
    }

    sResult = sPath;

    return sResult;
}
//---------------------------------------------------------------------------

// Converts the UNIX path to the Windows path
string get_windows_path(string sPath)
{
int i = 0;
string sResult = "";

    for(i = 0; i < int(sPath.length()); i++)
    {
        if(sPath[i] == '/')
        {
            sPath[i] = '\\';
        }
    }

    sResult = sPath;

    return sResult;
}
//---------------------------------------------------------------------------

// Completes the path if necessary
string complete_path(string sPath)
{
int iPathLength = 0;
bool bUNIX = true;
string sResult = "";

    iPathLength = sPath.length();
    if(iPathLength > 0)
    {
        if(int(sPath.find("\\")) >= 0 )
        {
            bUNIX = false;
        }

        if(bUNIX == true)
        {   // UNIX
            if(sPath[iPathLength-1] != '/')
            {
                sPath += "/";
            }
        }
        else
        {   // Windows
            if(sPath[iPathLength-1] != '\\')
            {
                sPath += "\\";
            }
        }
    }

    sResult = sPath;

    return sResult;
}
//---------------------------------------------------------------------------

// Set the current working directory
string set_working_dir(string sDestExe)
{
string sResult = "";
string sWorkingDir = "";
string sPath = "";
string sTempUnixPath = "";
char cWorkingDir[STRING_LENGTH] = {};

    //If something goes wrong while getting the working dir
    sResult = sDestExe;

    //Get current application directory
	GetModuleFileNameA(NULL, cWorkingDir, STRING_LENGTH);
    sWorkingDir = cWorkingDir;

    if(sDestExe.length() > 0 && sWorkingDir.length() > 0)
    {
        //Repair the path if needed
        sWorkingDir = repair_cygwin_path(sWorkingDir);

        //Convert Windows path to UNIX path
        sWorkingDir = get_unix_path(sWorkingDir);

        //Fix me:
        //Create filesystem::path object ( works in Cygwin only with UNIX path )
        filesystem::path p(sWorkingDir);

        //Get file path
        sPath = p.parent_path();

        //Completes the path if necessary
        sPath = complete_path(sPath);

        //Convert UNIX path to Windows path
        sPath = get_windows_path(sPath);

        //Build the filename
        sResult = "\"" + sPath + sDestExe + "\"";
    }

    return sResult;
}
//---------------------------------------------------------------------------

// Create the temp file
string create_temp_file(string sFilename)
{
string sTempPath = "";
string sFilenameWoPath = "";
string sTempUnixFilename = "";
string sTempFilename = "";

    //Fix me:
    //The result of getenv is not correct in cygwin environment.
    //
    //Example: getenv("TEMP");
    //
    //Result : /cygdrive/c/path/to/temp
    //But should be: c:\path\to\temp
    //
    //Get environment variable TEMP
    sTempPath = getenv("TEMP");

    //Repair the path if needed
    sTempPath = repair_cygwin_path(sTempPath);

    //Completes the path if necessary
    sTempPath = complete_path(sTempPath);

    //Convert Windows path to UNIX path
    sTempUnixFilename = get_unix_path(sFilename);

    //Fix me:
    //Create filesystem::path object ( works in Cygwin only with UNIX path )
    filesystem::path p(sTempUnixFilename);

    //Get filename without path
    sFilenameWoPath = p.filename();

    //Create filename
    sTempFilename = sTempPath + sFilenameWoPath;

    //Create temp file at destination
    if(CopyFile(sFilename.c_str(), sTempFilename.c_str(), false) == false)
    {
        die("Could not create temp file");
    }

    //Return file path
    return sTempFilename;
}
//---------------------------------------------------------------------------

// Delete file
bool delete_temp_file(string sFilename)
{
    return DeleteFile(sFilename.c_str());
}
//---------------------------------------------------------------------------

// Execute a command and get the results.
string ExecCmd(string sCMD /* [in] command to execute */ )
{
int i = 0;
string strResult = "";
char cmd[STRING_LENGTH] = {};
HANDLE hPipeRead = 0;
HANDLE hPipeWrite = 0;

    strcpy(cmd, sCMD.c_str());

    SECURITY_ATTRIBUTES saAttr = {sizeof(SECURITY_ATTRIBUTES)};
    saAttr.bInheritHandle = TRUE; // Pipe handles are inherited by child process.
    saAttr.lpSecurityDescriptor = NULL;

    // Create a pipe to get results from child's stdout.
    if (!CreatePipe(&hPipeRead, &hPipeWrite, &saAttr, 0))
        return strResult;

    STARTUPINFOA si = {sizeof(STARTUPINFOA)};
    si.dwFlags     = STARTF_USESHOWWINDOW | STARTF_USESTDHANDLES;
    si.hStdOutput  = hPipeWrite;
    si.hStdError   = hPipeWrite;
    si.wShowWindow = SW_HIDE; // Prevents cmd window from flashing.
                              // Requires STARTF_USESHOWWINDOW in dwFlags.

    PROCESS_INFORMATION pi = { 0 };

    BOOL fSuccess = CreateProcessA(NULL, (LPSTR)cmd, NULL, NULL, TRUE, CREATE_NEW_CONSOLE, NULL, NULL, &si, &pi);
    if (! fSuccess)
    {
        CloseHandle(hPipeWrite);
        CloseHandle(hPipeRead);
        return strResult;
    }

    bool bProcessEnded = false;
    for (; !bProcessEnded ;)
    {
        // Give some timeslice (50 ms), so we won't waste 100% CPU.
        bProcessEnded = WaitForSingleObject( pi.hProcess, 50) == WAIT_OBJECT_0;

        // Even if process exited - we continue reading, if
        // there is some data available over pipe.
        char buf[STRING_LENGTH] = {};
        DWORD dwRead = 0;
        DWORD dwAvail = 0;

        for (;;)
        {
            if (!PeekNamedPipe(hPipeRead, NULL, 0, NULL, &dwAvail, NULL))
                break;

            if (!dwAvail) // No data available, return
                break;

            if (!ReadFile(hPipeRead, buf, min( sizeof(buf) - 1, (long unsigned int )dwAvail ), &dwRead, NULL) || !dwRead)
                // Error, the child process might ended
                break;

            //Get results
            strResult += buf;

            //Clean up entire buffer to remove trash
            for(i = 0; i < STRING_LENGTH; i++){buf[i] = 0;}
        }
    }

    CloseHandle(hPipeWrite);
    CloseHandle(hPipeRead);
    CloseHandle(pi.hProcess);
    CloseHandle(pi.hThread);

    return strResult;
}
//---------------------------------------------------------------------------

//auto detect codec and bitrate from input stream, generate ffmpeg options for audio encoder
string audio_encode_options(string sFilename)
{
string sParam = "";
string sCMD = "";
string sResult = "";
string sValue = "";
string sCodecName = "";
string sProbeResult = "";
string sTempResult = "";
string sBitRate = "";
string sData = "";
int i = 0;
int iPos = 0;

    //ffprobe -v error -print_format compact -show_streams "$1"
    sParam = "-print_format compact -show_streams \"" + sFilename + "\"";
    sCMD = g_sFFProbe + " " + g_sFFMPEG_VERBOSE + " " + sParam;

    //Execute and get results
    sProbeResult = ExecCmd(sCMD);

    //Parse to find the audio line
    for(i = 0; i < int(sProbeResult.length()); i++)
    {
        if(sProbeResult[i] != 0x0d && sProbeResult[i] != 0x0a)
        {
            sTempResult += sProbeResult[i];
        }

        if(sProbeResult[i] == 0x0d)
        {
            if(int(sTempResult.find("codec_type=audio")) > 0)
            {
                sResult = sTempResult;
                break;
            }
            else
            {
                sTempResult = "";
            }
        }
    }

    //Parse to find codec_name and bit_rate
    if(sResult.length() > 0)
    {
        sValue = "codec_name=";
        iPos = int(sResult.find(sValue));
        if(iPos > 0)
        {
            iPos += sValue.length();
            while(sResult[iPos] != '|')
            {
                sCodecName += sResult[iPos];
                iPos++;
            }
        }

        sValue = "bit_rate=";
        iPos = int(sResult.find(sValue));
        if(iPos > 0)
        {
            iPos += sValue.length();
            while(sResult[iPos] != '|')
            {
                sBitRate += sResult[iPos];
                iPos++;
            }
        }

        // opus encoder is experimental, ffmpeg recommends libopus for encoding
        if (sCodecName == "opus")
        {
            sCodecName = "libopus";
        }

        sResult = "-c:a " + sCodecName;

        if (sBitRate != "N/A")
        {
            sResult += " -ab " + sBitRate;
        }
    }

    return sResult;
}
//---------------------------------------------------------------------------

// count number of audio and video streams, typical output: "audio=1:video=1"
string audio_video_stream_count(string sFilename)
{
string sParam = "";
string sCMD = "";
string sResult = "";
char cResult[STRING_LENGTH] = {};
int iVideoCodec = 0;
int iAudioCodec = 0;

    //ffprobe -v error -print_format compact -show_streams "$1"
    sParam = "-print_format compact -show_streams \"" + sFilename + "\"";
    sCMD = g_sFFProbe + " " + g_sFFMPEG_VERBOSE + " " + sParam;

    //Execute and get results
    sResult = ExecCmd(sCMD);

    if(sResult.length() > 0)
    {
        if(int(sResult.find("codec_type=audio")) > 0)
        {
            iAudioCodec = 1;
        }

        if(int(sResult.find("codec_type=video")) > 0)
        {
            iVideoCodec = 1;
        }

        sprintf(cResult, "audio=%d:video=%d", iAudioCodec, iVideoCodec);

        sResult = cResult;
    }

    return sResult;
}
//---------------------------------------------------------------------------

// get file extension
string extension(string sFilename)
{
string sResult = "";
string sTempUnixFilename = "";

    //Fix me:
    //filesystem::path in Cygwin only works with /
    sTempUnixFilename = get_unix_path(sFilename);

    filesystem::path p(sTempUnixFilename);

    sResult = p.extension();

    return sResult;
}
//---------------------------------------------------------------------------

// add watermark to file
void add_watermark(string sInFile, string sOutFile, string sHash, string sARGS)
{
string sCMD = "";
string sParam = "";
string sExtIn = "";
string sExtOut = "";
string sResult = "";
string sStreamCount = "";
string sEncodeOptions = "";
string sTempFilenameVideo = "";
string sTempFilenameAudio = "";
string sTempFilenameAudioWM = "";

    // check file extensions
    sExtIn = extension(sInFile);
    sExtOut = extension(sOutFile);
    if(sExtIn != sExtOut)
    {
        die("input/output extension must match ('" + sExtIn + "' vs. '" + sExtOut + "')");
    }

    // check audio/video stream count
    sStreamCount = audio_video_stream_count(sInFile);
    if(sStreamCount != "audio=1:video=1")
    {
        printf("videowmark: detected input file stream count: %s\n", sStreamCount.c_str());
        die ("input file must have one audio stream and one video stream");
    }

    // create tmpfiles
    sTempFilenameVideo = create_temp_file(sInFile);

    // create audio tmpfilename ( just the name, not the file )
    sTempFilenameAudio = sTempFilenameVideo + ".wav";

    // create wm_audio tmpfilename ( just the name, not the file )
    sTempFilenameAudioWM = sTempFilenameVideo + "_wm.wav";

    // get audio as wav
    //ffmpeg $FFMPEG_VERBOSE -y -i "$in_file" -f wav -rf64 always "$wav"
    sParam = "-y -i \"" + sTempFilenameVideo + "\" -f wav -rf64 always \"" + sTempFilenameAudio + "\"";
    sCMD = g_sFFMpeg + " " +  g_sFFMPEG_VERBOSE + " " + sParam;

    //Execute and get results
    sResult = ExecCmd(sCMD);

    //Show detailed FFMpeg results
    if(g_sFFMPEG_VERBOSE == "-v info"){printf("%s\n", sResult.c_str());}

    if(sResult.length() > 0)
    {
        if(int(sResult.find("Error")) > 0)
        {
            //Clean up
            delete_temp_file(sTempFilenameVideo);
            delete_temp_file(sTempFilenameAudio);
            delete_temp_file(sTempFilenameAudioWM);

            die("extracting audio from video failed (ffmpeg)");
        }
    }

    // add watermark
    // audiowmark add "${AUDIOWMARK_ARGS[@]}" "$orig_wav" "$wm_wav" "$bits" --set-input-label "$in_file" --set-output-label "$out_file" --output-format rf64"
    sEncodeOptions = audio_encode_options(sInFile);
    if(g_iQuiet == 0){printf("Audio Codec:  %s\n", sEncodeOptions.c_str());}
    sParam = "add " + sARGS + " \"" + sTempFilenameAudio + "\" \"" + sTempFilenameAudioWM + "\" " + sHash + " --set-input-label \"" + sInFile + "\" --set-output-label \"" + sOutFile + "\" --output-format rf64";
    sCMD = g_sAudiowmark + " " + sParam;

    //Execute and get results
    sResult = ExecCmd(sCMD);

    if(sResult.length() > 0)
    {
        if(int(sResult.find("Error")) > 0)
        {
            //Show detailed cause of error
            if(g_sFFMPEG_VERBOSE == "-v info"){printf("%s\n", sResult.c_str());}

            //Clean up
            delete_temp_file(sTempFilenameVideo);
            delete_temp_file(sTempFilenameAudio);
            delete_temp_file(sTempFilenameAudioWM);

            die("watermark generation failed (audiowmark)");
        }
    }

    //Show Audiowmark results
    if(g_iQuiet == 0)
    {
        printf("%s\n", sResult.c_str());
    }

    // rejoin
    // ffmpeg $FFMPEG_VERBOSE -y -i "$in_file" -i "$wm_wav" -c:v copy $(audio_encode_options "$in_file") -map 0:v:0 -map 1:a:0 "$out_file"
    sParam = "-y -i \"" + sTempFilenameVideo + "\" -i \"" + sTempFilenameAudioWM + "\" -c:v copy " + sEncodeOptions + " -map 0:v:0 -map 1:a:0 \"" + sOutFile + "\"";
    sCMD = g_sFFMpeg + " " +  g_sFFMPEG_VERBOSE + " " + sParam;

    //Execute and get results
    sResult = ExecCmd(sCMD);

    //Show detailed FFMpeg results
    if(g_sFFMPEG_VERBOSE == "-v info"){printf("%s\n", sResult.c_str());}

    if(sResult.length() > 0)
    {
        if(int(sResult.find("Error")) > 0)
        {
            //Clean up
            delete_temp_file(sTempFilenameVideo);
            delete_temp_file(sTempFilenameAudio);
            delete_temp_file(sTempFilenameAudioWM);

            die("merging video and watermarked audio failed (ffmpeg)");
        }
    }

    //Clean up
    delete_temp_file(sTempFilenameVideo);
    delete_temp_file(sTempFilenameAudio);
    delete_temp_file(sTempFilenameAudioWM);
}
//---------------------------------------------------------------------------

// get watermark from file
void get_watermark(string sFilename, string sARGS)
{
string sParam = "";
string sCMD = "";
string sTempPath = "";
string sResult = "";
string sTempFilenameVideo = "";
string sTempFilenameAudio = "";

    // check audio/video stream count
    sResult = audio_video_stream_count(sFilename);

    if(sResult.length() > 0)
    {
        if(sResult != "audio=1:video=1" )
        {
            printf("videowmark: detected input file stream count: %s\n", sResult.c_str());
            die("input file must have one audio stream and one video stream");
        }
    }

    // create video tmpfile
    sTempFilenameVideo = create_temp_file(sFilename);

    // create audio tmpfilename ( just the name, not the file )
    sTempFilenameAudio = sTempFilenameVideo + ".wav";

    // get audio as wav
    //ffmpeg $FFMPEG_VERBOSE -y -i "$in_file" -f wav -rf64 always "$wav" || die "extracting audio from video failed (ffmpeg)"
    sParam = "-y -i \"" + sTempFilenameVideo + "\" -f wav -rf64 always \"" + sTempFilenameAudio + "\"";
    sCMD = g_sFFMpeg + " " + g_sFFMPEG_VERBOSE + " " + sParam;

    //Execute and get results
    sResult = ExecCmd(sCMD);

    //Show detailed FFMpeg results
    if(g_sFFMPEG_VERBOSE == "-v info"){printf("%s\n", sResult.c_str());}

    if(sResult.length() > 0)
    {
        if(int(sResult.find("Error")) > 0)
        {
            //Clean up
            delete_temp_file(sTempFilenameVideo);
            delete_temp_file(sTempFilenameAudio);

            die("extracting audio from video failed (ffmpeg)");
        }
    }

    // get watermark
    //audiowmark get "${AUDIOWMARK_ARGS[@]}" "$wav" || die "retrieving watermark from audio failed (audiowmark)"
    sParam = "get " + sARGS + " \"" + sTempFilenameAudio + "\"";
    sCMD = g_sAudiowmark + " " + sParam;

    //Execute and get results
    sResult = ExecCmd(sCMD);

    if(sResult.length() > 0)
    {
        if(int(sResult.find("Error")) > 0)
        {
            //Show detailed cause of error
            if(g_sFFMPEG_VERBOSE == "-v info"){printf("%s\n", sResult.c_str());}

            //Clean up
            delete_temp_file(sTempFilenameVideo);
            delete_temp_file(sTempFilenameAudio);

            die("retrieving watermark from audio failed (audiowmark)");
        }
    }

    //Show Audiowmark results
    printf("%s\n", sResult.c_str());

    //Clean up
    delete_temp_file(sTempFilenameVideo);
    delete_temp_file(sTempFilenameAudio);
}
//---------------------------------------------------------------------------

void show_version_and_exit()
{
    printf("%s\n", g_sVersion.c_str());
    exit(0);
}
//---------------------------------------------------------------------------

void show_help_and_exit()
{
    printf(
           "usage: videowmark <command> [ <args>... ]\n"
           "\n"
           "Commands:\n"
           "  * create a watermarked video file with a message\n"
           "    videowmark add <input_video> <watermarked_video> <message_hex>\n"
           "\n"
           "  * retrieve message\n"
           "    videowmark get <watermarked_video>\n"
           "\n"
           "Global options:\n"
           "  --strength <s>            set watermark strength\n"
           "  --key <file>              load watermarking key from file\n"
           "  -q, --quiet               disable information messages\n"
           "  -v, --verbose             enable ffmpeg verbose output\n"
           "  --version                 show the current videowmark version\n"
           "  --detect-speed            detect and correct replay speed difference\n"
           "  --detect-speed-patient    slower, more accurate speed detection\n"
           );

    exit(0);
}
//---------------------------------------------------------------------------

int main(int argc , char *argv[])
{
string sResult = "";
string sAction = "";
string sKeyfile = "";
string sFilename = "";
string sInFile = "";
string sOutFile = "";
string sHash = "";
string sARGS = "";
int i = 0;

    // Set the the current working directory
    g_sFFProbe = set_working_dir(g_sFFProbe);
    g_sFFMpeg = set_working_dir(g_sFFMpeg);
    g_sAudiowmark = set_working_dir(g_sAudiowmark);

    if(argc > 1)
    {
        // Get all args
        for(i = 1; i < argc; i++)
        {
            if(string(argv[i]) == "-v" || string(argv[i]) == "--verbose")
            {
                g_sFFMPEG_VERBOSE = "-v info";
            }
            else if(string(argv[i]) == "--version")
            {
                show_version_and_exit();
            }
            else if(string(argv[i]) == "-q" || string(argv[i]) == "--quiet")
            {
                sARGS += " -q";
                g_iQuiet = 1;
            }
            else if(string(argv[i]) == "--detect-speed" || string(argv[i]) == "--detect-speed-patient")
            {
                sARGS += " " + string(argv[i]);
            }
            else if(string(argv[i]) == "-h" || string(argv[i]) == "--help")
            {
                show_help_and_exit();
            }
            else if(string(argv[i]) == "--key")
            {
                if(argc >= i+1 )
                {
                    sARGS += " " + string(argv[i]) + " " + string(argv[i+1]);
                    i++;
                }
                else
                {
                    die("videowmark: error parsing command line arguments (use videowmark -h)");
                }
            }
            else if(string(argv[i]) == "--strength")
            {
                if(argc >= i+1 )
                {
                    sARGS += " " + string(argv[i]) + " " + string(argv[i+1]);
                    i++;
                }
                else
                {
                    die("videowmark: error parsing command line arguments (use videowmark -h)");
                }
            }
            else if(string(argv[i]) == "add" || string(argv[i]) == "get" || string(argv[i]) == "probe")
            {
                sAction = string(argv[i]);
            }
            else if(sInFile.length() == 0 && string(argv[i]).length() > 0)
            {
                sInFile = string(argv[i]);
            }
            else if(sOutFile.length() == 0 && string(argv[i]).length() > 0)
            {
                sOutFile = string(argv[i]);
            }
            else if(sHash.length() == 0 && string(argv[i]).length() > 0)
            {
                sHash = string(argv[i]);
            }
            else
            {
                //
            }
        }

        // Get and execute action
        if(sAction == "add" && sInFile.length() > 0 && sOutFile.length() > 0 && sHash.length() > 0)
        {
            add_watermark(sInFile, sOutFile, sHash, sARGS);
        }
        else if(sAction == "get" && sInFile.length() > 0)
        {
            get_watermark(sInFile, sARGS);
        }
        else if(sAction == "probe" && sInFile.length() > 0)
        {
            printf("%s %s\n", sInFile.c_str(),audio_encode_options(sInFile).c_str());
        }
        else
        {
            printf("videowmark: error parsing command line arguments (use videowmark -h)\n");
        }
    }
    else
    {
        show_help_and_exit();
    }

    return 0;
}
//---------------------------------------------------------------------------
