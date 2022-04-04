import sys


with open(sys.argv[1]) as f:
	lines = f.readlines()
	time = 0
	count = 0
	is_true = 0
	for line in lines:
		if "time_taken" in line:
			a = line.split(' ')
			time += float(a[3])
			count += 1
			if a[1] == 'True':
				is_true += 1
	print('avg_time:', time/count, 'is_true:', is_true)

