echo "########### Setting up env for expresso ##############"

system=$(uname --all)
if [[ "$system" == *"WSL2"* ]]; then
   echo "In wsl2"
   conda activate py37_coffea_hep
elif [[ "$system" == *"lxslc"* ]]; then
    echo "In lxslc"
    conda activate expresso
elif [[ "$system" == *"lxplus"* ]]; then
    echo "In lxplus"
    source /cvmfs/sft.cern.ch/lcg/views/dev4cuda/latest/x86_64-centos7-gcc8-opt/setup.sh
    pip install objprint
    pip install pprintpp
else
    echo "Are you in supported node?"
fi      

chmod +x plot+.py
chmod +x expresso.py
pip install -e .
chmod +x modules/createJSON.py


ehelp () {
    python expresso.py --help
}


phelp () {
    python plot+.py --help
}

ana () {
    ls Analysis/
}

testana () {
    ./expresso.py --Samples Analysis/barebones/samples.txt --Analysis barebones --NumberOfTasks 2 --Debug --SaveRoot
}
testanatight () {
    ./expresso.py --Samples Analysis/barebones/samples.txt --Analysis barebones --NumberOfTasks 2 --Debug --SaveRoot --AnalysisPoint tight_ele_tight_mu
}
testplot () {
    python Analysis/barebones/plot.py
    }

lumical () {
    if [ "$1" == "-h" ] ; then
	echo "Usage: lumical <path to data json> <full or mask> "
	echo "Usage: lumical blabla/DoubleMuon.json mask"
    return
    fi

    rm -rf Analysis/LumiCal/lumi.txt; ./expresso.py --NumberOfTasks 1000000 --Analysis Analysis/LumiCal/ --PassOptions $2 --Sample $1 --ChunkSize 1000000;awk '{ sum += $1 } END { printf "Luminosity is:  %d pb⁻¹\n", sum}' Analysis/LumiCal/lumi.txt;rm -rf Analysis/LumiCal/lumi.txt
    }
