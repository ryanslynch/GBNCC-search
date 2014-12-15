# For a Moab batch system, this module submits scripts
# and checks to see when they finish. It will run under
# Linux and Windows (where we have Moab over CCS).

import unittest
import re, os, time
import platform
from subprocess import Popen, PIPE
from xml.dom.minidom import parseString


_submission_script = "qsub"
_debug = True

def submit_job(filename, options=""):
    '''Either you get a job id, and your job will run, or you don't
    and your job won't run, so you check out and err to see why.'''
    global _submission_script
    process = Popen(_submission_script+" "+options+" "+filename, shell=True,
                    stdout=PIPE, stderr=PIPE)
    returncode = process.returncode
    (out,err) = process.communicate()
    return (find_jobid(out), out, err)


def find_jobid(message):
    '''Finds the jobid in the output from qsub. The job id can either
    be just a number or Moab.number.'''
    # Optional \r because Windows python2.4 can't figure out \r\n is newline
    match=re.search("^((Moab.)?(\d+))[\r]?",message,re.IGNORECASE|re.MULTILINE)
    if match:
        return match.group(1)
    return None


def get_all_jobs():
    '''Returns all jobs in queue as dictionary of jobs by job id
    containing dictionary for each job as key-value.'''
    process = Popen('showq --format=XML', shell=True,
                    stdout=PIPE, stderr=PIPE)
    (out,err) = process.communicate()
    if not out: return None
    try:
        # There is an odd Windows python2.4 bug here
        noneAtEnd=re.search("None\s*$",out)
        if noneAtEnd:
            out = out[:noneAtEnd.start()]
        dom=parseString(out)
        jobs=dom.getElementsByTagName('job')
        all_job_props={}
        for job in jobs:
            job_props={}
            for prop in job.attributes.keys():
                job_props[prop]=job.attributes[prop].nodeValue
            all_job_props[job_props['JobID']]=job_props
    except:
        return None
    return all_job_props
    


def get_job_state(jobid):
    global _debug
    process = Popen('checkjob -v --format=XML %s' % (jobid,), shell=True,
                    stdout=PIPE, stderr=PIPE)
    (out,err) = process.communicate()
    if not out: return None
    try:
        # There is an odd Windows python2.4 bug here
        noneAtEnd=re.search("None\s*$",out)
        if noneAtEnd:
            out = out[:noneAtEnd.start()]
        dom=parseString(out)
        jobs=dom.getElementsByTagName('job')
        job=jobs[0]
        job_props={}
        for prop in job.attributes.keys():
            job_props[prop]=job.attributes[prop].nodeValue
    except:
        return None
    return job_props



def is_job_done(jobid):
    '''Assume job did start, so if the checkjob never heard of this job
    it must already have finished.'''
    state=get_job_state(jobid)
    if state:
        if state['State']=='Completed':
            return True
        else:
            return False
    return True



def are_jobs_done(jobids):
    jobs=get_all_jobs()
    if not jobs: return True
    done=True
    for jobid in jobids:
        if jobs.has_key(jobid):
            if jobs[jobid]['State']!='Completed':
                done=False
    return done


def cancel_job(jobid):
    process = Popen('mjobctl -c %s' % (jobid,), shell=True,
                    stdout=PIPE, stderr=PIPE)
    (out,err)=process.communicate()


class SubmitTestCase(unittest.TestCase):
    def testWrongExe(self):
        global _submission_script
        ss=_submission_script
        _submission_script = 'blah'
        self.assert_(submit_job('job.sh')[0]==None)
        _submission_script = ss
        
    def testFindJobId(self):
        self.assert_('12344'==find_jobid(os.linesep.join(
            ["blahblah","12344"])))
        self.assert_('12344'==find_jobid(os.linesep.join(
            ["blahblah","12344"])))
        self.assert_('Moab.12344'==find_jobid(os.linesep.join(
            ["blahblah","Moab.12344"])))
        self.assert_(None==find_jobid(""))
        self.assert_(None==find_jobid(os.linesep.join(
            ["Error my friend 2367","more error"])))

    def testGetJobStateFail(self):
        res=get_job_state("nosuchjob")
        self.assert_(res==None)


    def testIsDoneLongDone(self):
        self.assert_(is_job_done('123')==True)


    def testGetJobState(self):
        jobs=get_all_jobs()
        res=get_job_state(jobs.keys()[0])
        self.assert_(res!=None)


    def testGetAllJobs(self):
        jobs=get_all_jobs()
        self.assert_(None != jobs)
        for jobId in jobs.keys():
            pass #print jobId
            for key in jobs[jobId].keys():
                pass #print "\t",key,"\t",jobs[jobId][key]


    def testJobsNotDone(self):
        done = are_jobs_done('blah')
        self.assert_(done==True)


    def testJobsDone(self):
        done = are_jobs_done(['1003758',])
        self.assert_(done==True)


class ActuallyDoesSubmitTestCases(unittest.TestCase):
    def testSubmit(self):
        if platform.system()=='Windows': return
        script = open('z.sh',"w")
        print >>script, """
#PBS -A dal16_0001
#PBS -l nodes=1,walltime=10:00
#PBS -q v4dev
#PBS -N clans2win1"""
        script.close()
        (jobid,out,err)=submit_job('z.sh')
        self.assert_(jobid)
        os.unlink('z.sh')


    def testBadSubmitQueue(self):
        script = open('z.sh',"w")
        print >>script, """
#PBS -A dal16_0001
#PBS -l nodes=1,walltime=10:00
#PBS -q v5000
#PBS -N clans2win1"""
        script.close()
        (jobid,out,err)=submit_job('z.sh')
        if None==jobid:
            print 'out',out
            print 'err',err
        self.assert_(jobid==None)
        os.unlink('z.sh')


    def testWatchSubmit(self):
        if platform.system()=='Windows': return
        script = open('z.sh',"w")
        print >>script, """
#PBS -A dal16_0001
#PBS -l nodes=1,walltime=10:00
#PBS -q v4dev
#PBS -N clans2win1"""
        script.close()
        (jobid,out,err)=submit_job('z.sh')
        self.assert_(jobid)

        while not is_job_done(jobid):
            res=get_job_state(jobid)
            print res['State'], res['EState']
            time.sleep(5)
        
        os.unlink('z.sh')


    def testWatchSubmitWindows(self):
        script = open('z.sh',"w")
        print >>script, """
#PBS -A dal16_0001
#PBS -l nodes=1,walltime=10:00
#PBS -q NORMAL
#PBS -N clans2win1"""
        script.close()
        (jobid,out,err)=submit_job('z.sh')
        self.assert_(jobid)

        while not is_job_done(jobid):
            res=get_job_state(jobid)
            print res['State'], res['EState']
            time.sleep(5)
        
        os.unlink('z.sh')



def suite():
    s=unittest.TestSuite()
    s.addTest(unittest.TestLoader().loadTestsFromTestCase(SubmitTestCase))
    s.addTest(unittest.TestLoader().loadTestsFromTestCase(
        ActuallyDoesSubmitTestCases))
    return s

if __name__ == "__main__":
    unittest.TextTestRunner(verbosity=2).run(suite())
