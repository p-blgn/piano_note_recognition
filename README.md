# piano_note_recognition
A basic piano note recognition programme in C using Fourier transform as part of a final-year bachelor's course.

A function reads the bits of a WAV file and convert them into a table containing the values of each sample of the music file. Then a Fourier transform allows us to know the distribution in the frequency domain and hence detects the notes.

This programme uses needs GSL to work.
