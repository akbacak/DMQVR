
system('sh  ./Python/extractFrames.sh');

system('python Python/genNP.py');

system('python Python/genHashCodes.py');

system('python Python/genLabels.py');
