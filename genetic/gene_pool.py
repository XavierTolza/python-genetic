import numpy as np


class gene_pool(object):
    def __init__(self, alleles=None, indexes = None):
        self.__alleles = []
        self.__alleles_indexes = []
        if alleles is not None:
            if indexes is None:
                indexes = [None]*len(alleles)
            for allele,index in zip(alleles,indexes):
                self.add_allele(allele,index)

    @property
    def n_allelles(self):
        return len(self.__alleles)

    @property
    def alleles_indexes(self):
        return np.array(self.__alleles_indexes)

    @property
    def alleles(self):
        return np.array(self.__alleles)

    def add_allele(self,allele, index=None):
        self.__alleles.append(allele)

        if index is None:
            index = np.random.random()
        self.__alleles_indexes.append(index)

    def __gene_selector(self,start,stop):
        selector = self.alleles_indexes
        selector = np.logical_and(selector >= start,selector<=stop)
        return selector

    def __getslice__(self,i,j):
        selector = self.__gene_selector(i,j)
        return self.alleles[selector],self.__alleles_indexes[selector]

    def get_genes_outside_range(self,start,stop):
        selector = np.logical_not(self.__gene_selector(start, stop))
        return self.alleles[selector],self.alleles_indexes[selector]

    def copy(self):
        return self.__class__(self.alleles,self.alleles_indexes)

    def remove_slice(self,start,stop):
        selector = np.logical_not(self.__gene_selector(start, stop))
        self.__alleles = self.alleles[selector]
        self.__alleles_indexes = self.__alleles_indexes[selector]
