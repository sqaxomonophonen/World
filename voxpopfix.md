Building/running (on Linux):
 - Install SDL2 development libraries (`pkg-config --cflags sdl2` must succeed)
 - Build: `make voxpopfix`
 - Run: `./build/voxpopfix mm_left.wav` (WAV is included)
 - It takes a while to start the first time (see progress by running it in a terminal)

Controls:
 - Right click and drag: pan
 - Mouse wheel: zoom
 - Left click and drag: paint as broken
 - Shift + left click and drag: paint as not broken (or "erase broken")
 - Space: play "edited" audio from mouse cursor
 - D: play "original" audio from mouse cursor (as if nothing was marked broken, but audio is still synthesized)
 - R: render to disk
 - 1-6: toggle row visibility
 - Escape: quit

Notes / known problems:
 - Really slow software rendering; try changing the `DIVRES` define from 1 to 2 or higher in `test/voxpopfix.c` if you're on a high resolution display.
 - Crashes if you press space in the left out of data bounds region.
 - The editor is in `test/voxpopfix.c`, the rest is "library code" ("WORLD", SDL2, dr_wav)
 - Requires ~400MB of "analysis data" per minute of audio which is calculated once and stored on disk. It's likely that this could be compressed A LOT given that "WORLD" smells a bit like a voice compression research platform.
 - Besides usage of `mmap()` and friends, and GNU Make, the code ought to build/run on any platform that SDL2/dr_wav.h supports.
 
What is this software? Casey Muratori had some broken audio over at https://github.com/cmuratori/moriarty_audio and was asking for suggestions on how to fix it. My intuition was to transform the audio into a "smoother domain", so that missing/broken data could be replaced by linear interpolation of nearby "good" data (raw PCM data doesn't fit the "definition" of "smooth data" because it oscillates and cannot be lerped). I went with the domain that "WORLD" provides, which consists of 5ms frames containing f0 (fundamental frequency), spectrogram and aperiodicity data. It also means that the approach struggles a bit with polyphonic signals (e.g. 2+ simultaneous voices).

Rows from top to bottom:
 - "broken": what you have marked as broken.
 - "dry": not fully implemented, meant to crossfade back to original audio (when the synthesizer struggles)
 - "other spectrogram": settings chosen so that the problem is most visible
 - "f0 track": fundamental frequency (used in resynthesis)
 - "spectrogram track": rescaled (?) spectrogram (used in resynthesis)
 - "aperiodicity track": data used in resynthesis (don't really know what it does)
