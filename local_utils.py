'''
Created on 10 Dec 2014

@author: wbryant
'''

import sys
from time import time


def f_measure_tf(tp,tn,fp,fn):
    try:
        precision = tp / float(tp+fp)
        recall = tp / float(tp + fn)
        f_measure = 2*precision*recall/(precision + recall)
        accuracy = float(tp + tn) / (tp + tn + fp + fn)
        balanced_accuracy = 0.5 * (float(tp)/(tp + fn) + float(tn)/(tn+fp))  
        print("%12s = %1.3f" % ('Precision', precision))
        print("%12s = %1.3f" % ('Recall', recall))
        print("%12s = %1.3f" % ('F-measure', f_measure))
        print("%12s = %1.3f" % ('Accuracy', accuracy))
        print("%12s = %1.3f" % ('Bal. Accuracy', balanced_accuracy))
    except:
        f_measure = 0
    return f_measure


class ResultSet:
    """A container for a set of results that can return stats"""
    
    def __init__(self, tp, tn, fp, fn):
        
        self.tp = tp
        self.tn = tn
        self.fp = fp
        self.fn = fn
    
    def precision(self):
        try:
            return self.tp / float(self.tp+self.fp)
        except:
            return 0
    
    def recall(self):
        try:
            return self.tp / float(self.tp+self.fn)
        except:
            return 0
            
    def f_measure(self):
        precision = self.precision()
        recall = self.recall()
        try:
            return 2*precision*recall/(precision + recall)
        except:
            return 0
    
    def accuracy(self):
        try:
            return float(self.tp + self.tn) / (self.tp + self.tn + self.fp + self.fn)
        except:
            return 0
    
    def balanced_accuracy(self):
        try:
            return 0.5 * (float(self.tp)/(self.tp + self.fn) + float(self.tn)/(self.tn+self.fp))
        except:
            return 0
    
    def max_balanced_accuracy(self, num_tests_remaining):
        """Given a number of tests to be added to the results, what is the 
        maximum balanced accuracy that could be achieved?  
        """
        
        pass
    
    def calc_balanced_tfpn_addition(self, num_tests_remaining):
        denominator = self.fn + self.fp
        numerator = self.fn*self.tn - self.tp*self.fp + self.fn*num_tests_remaining
        
        x = round(numerator/float(denominator))
        y = num_tests_remaining - x
        
        return x, y
         
        
        
    def stats(self):
        if self.tp > 0:
            precision = self.precision()
            recall = self.recall()
            f_measure = self.f_measure()
            accuracy = self.accuracy()
            balanced_accuracy = self.balanced_accuracy()
            
            print("%12s = %1.3f" % ('Precision', precision))
            print("%12s = %1.3f" % ('Recall', recall))
            print("%12s = %1.3f" % ('F-measure', f_measure))
            print("%12s = %1.3f" % ('Accuracy', accuracy))
            print("%12s = %1.3f" % ('Bal. Accuracy', balanced_accuracy))
        else:
            print("There are no true positives.")



    
    

class loop_counter:
    """Use to track progress of a loop of known length."""
    
    def __init__(self, length, message = 'Entering loop', timed = False):
        self.stopped = False
        self.length = length
        self.num_done = 0
        self.next_progress = 1
        self.timed = timed
        if self.timed:
            self.time_0 = time()
        print("{}:".format(message))
        sys.stdout.write("\r - %d %%" % self.num_done)
        sys.stdout.flush()
    
    def step(self, number = None):
        self.num_done += 1
        if self.num_done > self.length:
            self.stop()
        if ((100 * self.num_done) / self.length) > self.next_progress:
            number = self.num_done-1
            if number is not None:
                if self.length > 100:
                    sys.stdout.write("\r - %d %%  (%d / %d)" % (self.next_progress, number, self.length))
                else:
                    percent_done = int(100*self.num_done/float(self.length))
                    sys.stdout.write("\r - %d %%  (%d / %d)" % (percent_done, number, self.length))
            else:
                sys.stdout.write("\r - %d %%" % self.next_progress)    
            if self.timed:
                time_n = time()
                sys.stdout.write(" - {} seconds".format(time_n-self.time_0))
            sys.stdout.flush()
            self.next_progress += 1
                
    def stop(self):
        if self.stopped == False:
            sys.stdout.write("\r - 100 %\n")
            self.stopped = True
            if self.timed:
                time_n = time()
                sys.stdout.write(" taking {} seconds.\n".format(time_n-self.time_0))
            sys.stdout.flush()
#         elif self.timed:
#             time_n = time()
#             sys.stdout.write(" - {} seconds overall ...".format(time_n-self.time_0))