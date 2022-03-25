import os
import sys

def prodigal_score_rbs(seq):
	s = seq[::-1]
	score = 0
	if 'ggagga' in (s[5:11],s[6:12],s[7:13],s[8:14],s[9:15],s[10:16]):
		score = 27
	elif 'ggagga' in (s[3:9],s[4:10]):
		score = 26
	elif 'ggagga' in (s[11:17],s[12:18]):
		score = 25
	elif 'ggagg' in (s[5:10],s[6:11],s[7:12],s[8:13],s[9:14],s[10:15]):
		score = 24
	elif 'ggagg' in (s[3:8],s[4:9]):
		score = 23
	elif 'gagga' in (s[5:10],s[6:11],s[7:12],s[8:13],s[9:14],s[10:15]):
		score = 22
	elif 'gagga' in (s[3:8],s[4:9]):
		score = 21
	elif 'gagga' in (s[11:16],s[12:17]) or 'ggagg' in (s[11:16],s[12:17]):
		score = 20
	elif 'ggacga' in (s[5:11],s[6:12],s[7:13],s[8:14],s[9:15],s[10:16]):
		score = 19
	elif 'ggatga' in (s[5:11],s[6:12],s[7:13],s[8:14],s[9:15],s[10:16]):
		score = 19
	elif 'ggaaga' in (s[5:11],s[6:12],s[7:13],s[8:14],s[9:15],s[10:16]):
		score = 19
	elif 'ggcgga' in (s[5:11],s[6:12],s[7:13],s[8:14],s[9:15],s[10:16]):
		score = 19
	elif 'ggggga' in (s[5:11],s[6:12],s[7:13],s[8:14],s[9:15],s[10:16]):
		score = 19
	elif 'ggtgga' in (s[5:11],s[6:12],s[7:13],s[8:14],s[9:15],s[10:16]):
		score = 19
	elif 'ggaaga' in (s[3:9],s[4:10]) or 'ggatga' in (s[3:9],s[4:10]) or 'ggacga' in (s[3:9],s[4:10]):
		score = 18
	elif 'ggtgga' in (s[3:9],s[4:10]) or 'ggggga' in (s[3:9],s[4:10]) or 'ggcgga' in (s[3:9],s[4:10]):
		score = 18
	elif 'ggaaga' in (s[11:17],s[12:18]) or 'ggatga' in (s[11:17],s[12:18]) or 'ggacga' in (s[11:17],s[12:18]):
		score = 17
	elif 'ggtgga' in (s[11:17],s[12:18]) or 'ggggga' in (s[11:17],s[12:18]) or 'ggcgga' in (s[11:17],s[12:18]):
		score = 17
	elif 'ggag' in (s[5:9],s[6:10],s[7:11],s[8:12],s[9:13],s[10:14]):
		score = 16
	elif 'gagg' in (s[5:9],s[6:10],s[7:11],s[8:12],s[9:13],s[10:14]):
		score = 16
	elif 'agga' in (s[5:9],s[6:10],s[7:11],s[8:12],s[9:13],s[10:14]):
		score = 15
	elif 'ggtgg' in (s[5:10],s[6:11],s[7:12],s[8:13],s[9:14],s[10:15]):
		score = 14
	elif 'ggggg' in (s[5:10],s[6:11],s[7:12],s[8:13],s[9:14],s[10:15]):
		score = 14
	elif 'ggcgg' in (s[5:10],s[6:11],s[7:12],s[8:13],s[9:14],s[10:15]):
		score = 14
	elif 'agg' in (s[5:8],s[6:9],s[7:10],s[8:11],s[9:12],s[10:13]):
		score = 13
	elif 'gag' in (s[5:8],s[6:9],s[7:10],s[8:11],s[9:12],s[10:13]):
		score = 13
	elif 'gga' in (s[5:8],s[6:9],s[7:10],s[8:11],s[9:12],s[10:13]):
		score = 13
	elif 'agga' in (s[11:15],s[12:16]) or 'gagg' in (s[11:15],s[12:16]) or 'ggag' in (s[11:15],s[12:16]):
		score = 12
	elif 'agga' in (s[3:7],s[4:8]) or 'gagg' in (s[3:7],s[4:8]) or 'ggag' in (s[3:7],s[4:8]):
		score = 11
	elif 'gagga' in (s[13:18],s[14:19],s[15:20]) or 'ggagg' in (s[13:18],s[14:19],s[15:20]) or 'ggagga' in (s[13:19],s[14:20],s[15:21]):
		score = 10
	elif 'gaaga' in (s[5:10],s[6:11],s[7:12],s[8:13],s[9:14],s[10:15]):
		score = 9
	elif 'gatga' in (s[5:10],s[6:11],s[7:12],s[8:13],s[9:14],s[10:15]):
		score = 9
	elif 'gacga' in (s[5:10],s[6:11],s[7:12],s[8:13],s[9:14],s[10:15]):
		score = 9
	elif 'ggtgg' in (s[3:8],s[4:9]) or 'ggggg' in (s[3:8],s[4:9]) or 'ggcgg' in (s[3:8],s[4:9]):
		score = 8
	elif 'ggtgg' in (s[11:16],s[12:17]) or 'ggggg' in (s[11:16],s[12:17]) or 'ggcgg' in (s[11:16],s[12:17]):
		score = 7
	elif 'agg' in (s[11:14],s[12:15]) or 'gag' in (s[11:14],s[12:15]) or 'gga' in (s[11:14],s[12:15]):
		score = 6
	elif 'gaaga' in (s[3:8],s[4:9]) or 'gatga' in (s[3:8],s[4:9]) or 'gacga' in (s[3:8],s[4:9]):
		score = 5
	elif 'gaaga' in (s[11:16],s[12:17]) or 'gatga' in (s[11:16],s[12:17]) or 'gacga' in (s[11:16],s[12:17]):
		score = 4
	elif 'agga' in (s[13:17],s[14:18],s[15:19]) or 'gagg' in (s[13:17],s[14:18],s[15:19]) or 'ggag' in (s[13:17],s[14:18],s[15:19]):
		score = 3
	elif 'agg' in (s[13:16],s[14:17],s[15:18]) or 'gag' in (s[13:16],s[14:17],s[15:18]) or 'gga' in (s[13:16],s[14:17],s[15:18]):
		score = 2
	elif 'ggaaga' in (s[13:19],s[14:20],s[15:21]) or 'ggatga' in (s[13:19],s[14:20],s[15:21]) or 'ggacga' in (s[13:19],s[14:20],s[15:21]):
		score = 2
	elif 'ggtgg' in (s[13:18],s[14:19],s[15:20]) or 'ggggg' in (s[13:18],s[14:19],s[15:20]) or 'ggcgg' in (s[13:18],s[14:19],s[15:20]):
		score = 2
	elif 'agg' in (s[3:6],s[4:7]) or 'gag' in (s[3:6],s[4:7]) or 'gga' in (s[3:6],s[4:7]):
		score = 1
	return score


