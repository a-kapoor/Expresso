#python Expresso.py --Sample SAMPLES/background_samples/central_UL/UL16_DY50.json --OutputName results --OutputFolder Analysis/chflip/output --ChunkSize 1000 --NumberOfTasks 4 --Analysis chflip -pre Analysis/chflip/preprocessor.py -plotter Analysis/chflip/plot.py --Xrootd ''

python Expresso.py --Sample SAMPLES/test.json --OutputName results --OutputFolder Analysis/chflip/output --ChunkSize 1000 --NumberOfTasks 4 --Analysis chflip -pre Analysis/chflip/preprocessor.py -plotter Analysis/chflip/plot.py --Xrootd '' --SaveRoot
