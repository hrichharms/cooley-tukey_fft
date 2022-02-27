from scipy.fftpack import fft
from math import e, pi, cos


if __name__ == "__main__":

	i = complex(0, 1)

	from sys import argv
	N = int(argv[1])

	t = []
	t.append(fft([0] * N))
	t.append(fft([1] * N))
	t.append(fft([1 - k % 2 * 2 for k in range(N)]))
	t.append(fft([e ** (16 * pi * i * k / N) for k in range(N)]))
	t.append(fft([cos(16 * pi * k / N) for k in range(N)]))
	t.append(fft([e ** (43/7 * i * 2 * pi * k / N) for k in range(N)]))
	t.append(fft([cos(43/7 * 2 * pi * k / N) for k in range(N)]))

	for k, results in enumerate(t):
		print(f"TEST {k + 1}:")
		print("|" + "|".join([f"({float(num.real)},{float(num.imag)})" for num in results]) + "|")
