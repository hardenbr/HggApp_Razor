#!/usr/bin/env python

import os
import sys
import getopt
import getpass

#queueString = "all.q@compute-2-4.local,all.q@compute-3-2.local,all.q@compute-3-3.local,all.q@compute-3-4.local,all.q@compute-3-7.local,all.q@compute-3-8.local,all.q@compute-3-9.local,all.q@compute-3-10.local,all.q@compute-3-11.local,all.q@compute-3-12.local"
queueString = "all.q@compute-3-2.local,all.q@compute-3-3.local,all.q@compute-3-4.local,all.q@compute-3-9.local,all.q@compute-3-10.local,all.q@compute-3-12.local"
#queueString = "all.q@compute-2-4.local,all.q@compute-3-5.local,all.q@compute-3-6.local,all.q@compute-3-7.local,all.q@compute-3-8.local,all.q@compute-3-11.local"
#queueString = "all.q@compute-2-4.local,all.q@compute-3-2.local,all.q@compute-3-7.local,all.q@compute-3-8.local"



def usage():
    print sys.argv[0]+" prepares jobs for the batch system and submits them"
    print sys.argv[0]+" [options] FileList "
    print "FileList is a text file containing the exact root path to each file on a new line"
    print "Options:"
    print "--queue=<queue>                   submits to queue <queue> (default: <8hn>)"
    print "--FilesPerJob=<#>                 submits <#> files per job (default: 10)"
    print "--DirName=<dir>                   puts job files in <dir> (default: ./batch)"
    print "--CMSSWVersion=<V>                specify the CMSSW Version to use, must specify the full name with CMSSW (default: CMSSW_5_2_3)"
    print "-6                                use CMSSW_6_2_0"
    print "--Arch=<A>                        specify the architecture for CMSSW (default:slc5_amd64_gcc462) "
    print "--submit                          submits the jobs **make sure you have tested first**"
    print "--noClear                         don't clear out the batch directory before making new jobs"
    print "--toy=N                           specify this if one is generating toys (don't need FileList, just create N jobs)"
    print "--MakeTmp=                        build the application in a temporary area and run from there"
    print ">>>> All of the following should be specified the same number of times!"
    print "--CMD=<opt>                       The exact command (with options and arguments) to execute (This option must be specied"
    print "                                      and can be specified multiple times for sequential commands) USE #IFL# TO SPECIFY THE INPUT LIST"
    print "                                      use #OF<num># (e.g. #OF0#,#OF11#,etc.) to specify the output files"
    print "--outputName=<name>               the base name of the output file (default: output), the number of the job will be inserted"
    print "                                      (defaults will be output_c1, output_c2, output_c3, etc.)"
    print "                                      only files specified here will be copied out of the /tmp area"
    print "--OutputDir=<dir>                 sends job outputs to <dir> (default: <jobdir>/out/)"

def insertString(pos,basestr,insstr):
    if pos <=0:
        return insstr+basestr
    elif pos >=len(basestr):
        return basestr+insstr
    else:
        return basestr[:pos]+insstr+basestr[pos:]

def main():
    try:
        opt, args = getopt.getopt(sys.argv[1:],"6",["queue=","FilesPerJob=","DirName=","OutputDir=",\
                                                    "submit","CMD=","outputName=","CMSSWVersion=","Arch=","FHADOOP","noClear","toy=",\
                                                    "MakeTmp"])
    except getopt.GetoptError, err:
            print str(err)
            print "Error"
            usage()
            sys.exit(2)

    #get required arguments

    wd = os.getcwd()

    #get options
    queue="8nh"
    FilesPerJob=10
    DirName=wd+"/batch"

    submitMode = False

    CMSSWVersion = "CMSSW_6_2_0"
    Arch = "slc5_amd64_gcc462"

    cmdStrings = []
    outputNames = []
    OutputDir=[]

    fHadoopMode = False
    noClear = False
    makeTmp = False
    nToy=-1
    for o,a in opt:
        if o=="--queue":
            queue = a
        elif o=="--FilesPerJob":
            FilesPerJob = int(a)
        elif o=="--DirName":
            DirName = a
        elif o == "--OutputDir":
            OutputDir.append(a)
        elif o == "--submit":
            submitMode = True
        elif o=="--CMD":
            print "Appending Command: "+a
            cmdStrings.append(a)
        elif o=="--outputName":
            outputNames.append(a)            
        elif o=="--CMSSWVersion":
            CMSSWVersion = a
        elif o=="--Arch":
            Arch = a
        elif o=="--FHADOOP":
            fHadoopMode=True
        elif o=="--noClear":
            noClear=True
        elif o=="--toy":
            nToy = int(a)
        elif o=="-6":
            CMSSWVersion = "CMSSW_6_2_0"
            Arch="slc5_amd64_gcc472"
        elif o=="--MakeTmp":
            makeTmp = True
        else:
            print "Invalid Option: %s with argument %s" % (o,a,)
            usage()
            sys.exit(1)

    if len(args)<1 and nToy < 0:
        usage()
        sys.exit(0)

    FileList = ""
    if nToy<0: FileList = args[0]

    if len(cmdStrings) == 0:
        print "Command String (--CMD) must be specified!\n"
        usage()
        sys.exit(0)

    while len(outputNames) < len(cmdStrings):
        outputNames.append("output_c%d.root" % (len(outputNames)+1,))

    while len(OutputDir) < len(cmdStrings):
        OutputDir.append("%s/out"%DirName)

    fNum=0
    if nToy<0:
        try:
            print "Opening: "+FileList
            FListFile = open(FileList)
            FList = FListFile.readlines()
        except:
            print "Error: could not open FileList \n"
            usage()
            raise
        FListFile.close()
    
        if len(FList)==0:
            print "No files specified \n"
            sys.exit(1)
        fNum=len(FList)

    HomeDir = os.environ['HOME']
    
    if not os.path.exists( "%s/%s"%(HomeDir,CMSSWVersion,) ):
        print "\n\nFATAL ERROR: you must setup %s/%s before submitting\n" % (HomeDir,CMSSWVersion,)
        sys.exit(0)

    ##setup the working dir
    if not noClear: os.system("rm -r %s" % DirName)
    os.system("mkdir -p %s" % DirName)
    os.system("mkdir -p %s/src" % DirName)
    os.system("mkdir -p %s/input" % DirName)
    os.system("mkdir -p %s/log" % DirName)
    for dir in OutputDir:
        if dir.startswith('/castor/cern.ch'):
            os.system("rfmkdir -p %s" % dir)
        elif dir.startswith('/eos'):
            os.system("cmsMkdir -p %s" % dir)
        elif dir.startswith('/mnt/hadoop') and fHadoopMode:
            os.system("hadoop fs -mkdir /%s" % dir.lstrip("/mnt/hadoop"))
        else:
            os.system("mkdir -p %s" % dir)

    NLists = 0

    if nToy <0:
        print "%d Files in input list" % (fNum,)
        while FilesPerJob*NLists < fNum:
            list = open("%s/input/input_%d.list" % (DirName,NLists,),'w')
            for i in range(FilesPerJob):
                try:  ## number of files may not be evenly divisible by FilesPerJob
                    list.write(FList[FilesPerJob*NLists+i].rstrip('\n')+"\n")  
                except:
                    pass
            list.close()
            NLists+=1

    else: ##toy mode
        NLists = nToy

    print "Will create %d jobs" % NLists

    for i in range(NLists):
        thisInput  = "/home/jhardenbrook/2013/RAZOR_DIPHOTON/HggApp_Razor/%s/input/input_%d.list" % (DirName,i,)            
        if nToy >=0: thisInput = ""
        script = []
        rmCommands=[]
        script.append("#!/bin/sh")
        script.append("#$ -S /bin/sh") # need this at caltech
        script.append("source /cvmfs/cms.cern.ch/cmsset_default.sh") # need this at caltech
#        script.append("cd %s/%s/src   # setup CMSSW" % (HomeDir,CMSSWVersion,) )

        script.append("cd /home/jhardenbrook/CMSSW_6_2_0")
        script.append("export HADOOP_CONF_DIR=/etc/hadoop")
        script.append("export SCRAM_ARCH=%s" % Arch)
        script.append("eval `scramv1 runtime -sh`")
        script.append("cd /wntmp/scratch/")
        script.append("rm -rf git_clone_%i" % i)
        script.append("mkdir git_clone_%i" % i)
        script.append("cd git_clone_%i " % i)
        script.append("git clone -b fake_sample https://github.com/hardenbr/HggApp_Razor.git")
        script.append("cd HggApp_Razor")
        script.append("cd /home/jhardenbrook/CMSSW_5_2_3")
        script.append("eval `scramv1 runtime -sh`")
        script.append("cd -")
        script.append("csh setup_fastjet_x64.sh")
        script.append("cd /home/jhardenbrook/CMSSW_6_2_0")
        script.append("eval `scramv1 runtime -sh`")
        script.append("cd -")
        script.append("make")
        script.append("cd /home/jhardenbrook/CMSSW_5_2_3")
        script.append("eval `scramv1 runtime -sh`")
        script.append("cd -")

        if makeTmp:
            script.append("make tmp")
        tmpList=""
        if fHadoopMode:
            files = open(thisInput).readlines()
            #tmpList = open("/wntmp/scratch/input_%d.list" % i,'w')
            script.append("rm /wntmp/scratch/input_%d.list" % i)
            for line in files:
                script.append("hadoop fs -get /%s /wntmp/scratch/" % line.lstrip('/mnt/hadoop').rstrip('\n'))
                #tmpList.write("/wntmp/scratch/%s\n" % os.path.basename(line).rstrip('\n') )
                tmpList = tmpList+"/wntmp/scratch/%s " % os.path.basename(line).rstrip('\n')
                script.append("echo /wntmp/scratch/%s >> /wntmp/scratch/input_%d.list" % (os.path.basename(line).rstrip('\n'),i,))
            #tmpList.close()
            thisInput = "/wntmp/scratch/input_%d.list" % i
            
        for cmd,dir,file in zip(cmdStrings,OutputDir,outputNames):
            #REPLACE THE SELECTOR RUNNING
            cmd = cmd.replace("#DIR#", "git_clone_%i/HggApp_Razor" % i)
            cmd = cmd.replace("#IFL#",thisInput)
            for j in range(len(outputNames)):
                cmd = cmd.replace("#OF%d#"%j, "/wntmp/scratch/"+insertString(outputNames[j].rfind('.'),outputNames[j],"_%d"%i))
            script.append(cmd)
            copyCMD = "cp"
            if dir.startswith('/castor/cern.ch'):
                copyCMD = "rfcp"
            elif dir.startswith('/eos'):
                copyCMD = "cmsStage"
            elif dir.startswith('/mnt/hadoop') and fHadoopMode:
                copyCMD = "hadoop fs -put"

            
                
#        for file in outputNames:
            if file != "NONE":
                if dir.startswith('/mnt/hadoop') and fHadoopMode:
                    script.append("hadoop fs -rm /%s/%s" % (dir.lstrip('/mnt/hadoop'),insertString(file.rfind('.'),file,"_%d"%i),) )
                    copy = "%s /wntmp/scratch/%s /%s" % (copyCMD,insertString(file.rfind('.'),file,"_%d"%i),dir.lstrip('/mnt/hadoop'),)
                else:
                    copy = "%s /wntmp/scratch/%s %s" % (copyCMD,insertString(file.rfind('.'),file,"_%d"%i),dir,)
                script.append(copy)
                rm = "rm /wntmp/scratch/%s" % insertString(file.rfind('.'),file,"_%d"%i)
                rmCommands.append(rm)
                rmCommands.append("rm -rf ../git_clone_%i" % i)
            
        if fHadoopMode:
            script.append("rm /%s" % tmpList)
            script.append("rm %s"% thisInput)


        script+=rmCommands
        scriptFile = "%s/src/submit_%d.src" % (DirName,i,)
        outScript = open(scriptFile,'w')
        for line in script:
            outScript.write(line.rstrip('\n')+'\n')
        outScript.close()

        if submitMode: # do the actual submission
            #queue = "all.q@compute-2-4.local,all.q@compute-3-2.local,all.q@compute-3-3.local,all.q@compute-3-4.local,all.q@compute-3-5.local,all.q@compute-3-6.local,all.q@compute-3-9.local,all.q@compute-3-10.local,all.q@compute-3-11.local,all.q@compute-3-12.local"
            queue = queueString
            #queue = "all.q@compute-2-2.local,all.q@compute-2-4.local,all.q@compute-3-2.local"
            logarea = wd + "/" + DirName
            bsub = "qsub -q  %s %s" % (queue,scriptFile)
            print bsub
            #os.system(bsub)
            
if __name__=="__main__":
    main()

