import numpy as np

class DrawPoint(object):
    def __init__(self,_i, _j, _bm):
        self.loc_i = _i
        self.loc_j = _j
        self.bm = _bm
        self.stream = []
        
        
    def set_influence(self,_sizex, _sizey, _depth):
        self.sizex = _sizex
        self.sizey = _sizey
        self.depth = _depth

        self.compute_influence()

    def set_influence(self,blocks):
        self.stream = blocks
        self.stream_len = len(self.stream)

    def compute_influence(self):
        self.stream = []
        self.level_stream = []
        for k in xrange(self.depth,self.bm.nz):
            level_stream = []
            for j in xrange((self.loc_j-self.sizey),(self.loc_j+self.sizey) + 1):
                if (j>=0 and j < self.bm.ny):
                    for i in xrange((self.loc_i-self.sizex),(self.loc_i+self.sizex) + 1):
                        if (i>=0 and i < self.bm.nx):
                            block_id = self.bm.get_blockid(i,j,k)
                            #cout << "adding block: " << i << ":" << j << ":" << k << ": blockid=" << block_id << endl;
                            level_stream += [block_id]

            self.level_stream += [level_stream]
            self.stream += level_stream
                            
        #check stream does not have duplicates
        self.stream_len = len(self.stream)
        
        ss = set(self.stream)
        assert self.stream_len == len(ss)

    def extraction(self,previousExtraction,units):
        
        extractionStart = previousExtraction
        extractionStop = extractionStart + units

        if (extractionStop > self.stream_len):
            extractionStop = self.stream_len

        # cout << "extractionStart: " << extractionStart << endl;
        # cout << "extractionStop: " << extractionStop << endl;


        n = extractionStop - extractionStart

        return list(self.stream[extractionStart:extractionStop])
