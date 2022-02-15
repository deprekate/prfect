import numpy as np

def find_median(sorted_list):
    indices = []
    list_size = len(sorted_list)
    median = 0
    if list_size % 2 == 0:
        indices.append(int(list_size / 2) - 1)  # -1 because index starts from 0
        indices.append(int(list_size / 2))
        median = (sorted_list[indices[0]] + sorted_list[indices[1]]) / 2
        pass
    else:
        indices.append(int(list_size / 2))
        median = sorted_list[indices[0]]
        pass
    return median, indices
    pass

def has_outlier(unsorted_list):
	if min(unsorted_list) < (np.mean(unsorted_list) - 2*np.std(unsorted_list)):
		return True
	else:
		return False
	sorted_list = sorted(unsorted_list)
	median, median_indices = find_median(sorted_list)
	q25, q25_indices = find_median(sorted_list[:median_indices[0]])
	q75, q75_indices = find_median(sorted_list[median_indices[-1] + 1:])
	iqr = q75 - q25
	lower = q25 - (iqr * 1.5)
	if sorted_list[0] < lower:
		return True
	else:
		return False
