# wavimg

A CLI tool for converting WAV audio files to PNG spectrogram images and reconstructing WAV files from PNG spectrograms.

## Features

* Generate a visual spectrogram (`.png`) from an audio file (`.wav`)
* Reconstruct an audio file (`.wav`) from its spectrogram image (`.png`)
* Adjust parameters for FFT - size, window length, hop length, decibel range, and Griffin-Lim iterations

## Usage

```bash
wavimg [OPTIONS] <input> <output>
```

| Option               | Description                         | Default |
| -------------------- | ----------------------------------- | ------- |
| `-c`, `--encode`     | Convert WAV to PNG spectrogram      |         |
| `-d`, `--decode`     | Convert PNG spectrogram to WAV      |         |
| `--n_fft <int>`      | FFT size                            | `2048`  |
| `--hop_length <int>` | Hop length                          | `512`   |
| `--win_length <int>` | Window length (0 = same as `n_fft`) | `0`     |
| `--min_db <float>`   | Minimum dB value                    | `-80.0` |
| `--max_db <float>`   | Maximum dB value                    | `0.0`   |
| `--gi <int>`         | Number of Griffin-Lim iterations    | `60`    |
| `--sr <int>`         | Sample rate (Hz)                    | `44100` |
| `-h`, `--help`       | Show help message and exit          |         |

## Examples

```bash
# Generate a PNG spectrogram from a WAV file
wavimg -c input.wav spectrogram.png

# Reconstruct a WAV file from a PNG spectrogram
wavimg -d spectrogram.png reconstructed.wav

# Generate a PNG spectrogram with custom parameters
wavimg -c input.wav spectrogram.png --n_fft=1024 --hop_length=256 --win_length=512
```

## Author

masuo / izawartka