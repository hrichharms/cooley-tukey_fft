from scipy.fftpack import fft
from math import e, pi, cos


def rf(x):
	return x if int(x) != x else int(x)


if __name__ == "__main__":

	i = complex(0, 1)

	from sys import argv
	N = int(argv[1])

	o = []
	t = []
	o.append([0.0] * N)
	t.append(fft(o[-1]))
	o.append([1.0] * N)
	t.append(fft(o[-1]))
	o.append([float(1 - k % 2 * 2) for k in range(N)])
	t.append(fft(o[-1]))
	o.append([e ** (16 * pi * i * k / N) for k in range(N)])
	t.append(fft(o[-1]))
	o.append([cos(16 * pi * k / N) for k in range(N)])
	t.append(fft(o[-1]))
	o.append([e ** (43/7 * i * 2 * pi * k / N) for k in range(N)])
	t.append(fft(o[-1]))
	o.append([cos(43/7 * 2 * pi * k / N) for k in range(N)])
	t.append(fft(o[-1]))

	for k, (original, results) in enumerate(zip(o, t)):
		print(f"TEST {k + 1}:")
		if type(original[0]) == float:
			print("\tINPUT DATA:\t[" + ",".join([f"({rf(num)},0)" for num in original]) + "]")
		else:
			print("\tINPUT DATA:\t[" + ",".join([f"({rf(num.real)},{rf(num.imag)})" for num in original]) + "]")
		print("\tOUTPUT DATA:\t[" + ",".join([f"({rf(num.real)},{rf(num.imag)})" for num in results]) + "]\n")
