# -*- coding: utf-8 -*-

#from sdfconf import sdf, mol2
try:
    import sdf, mol2
except ImportError:
    from sdfconf import sdf, mol2
from collections import OrderedDict as OrDi
import bisect as bi
import numpy

class Findable(object):
    
    def __init__(self, iterable=None, **kwargs):
        self._maindict = None #ordered dict jonka alkiot sisältävät koordinaatit listana, indeksi avaimena
        self._orders = [] #Lista joka tulee sisältämää dimensioiden määrän listapareja (tuplessa), jotka sisältävät kunkin koordinaattipisteen indeksin ja yhden dimension koordinaatin. Listat järjestetty kunkin dimension mukaan kasvavaan järjestykseen.
        self._thelen = 0 #Dimensioiden lukumäärä
        if iterable:
            self.set_iterable(iterable, kwargs.get('ignores', ['H'])) #Luodaan datarakenne
                
        else:
            self.set_iterable([]) #Tyhjä rakenne
            
    def set_iterable(self, iterable, ignores=[]):
        if type(iterable) in (list, tuple):
            pass
        elif type(iterable) in (sdf.Sdfmole,mol2.Mol2Mol):
            iterable = [list(atom[1]) for atom in iterable.atomsGenerator(ignores=ignores)]
        else:
            raise TypeError('Input must be of type list, tuple, Sdfmole, Mol2Mol.')
        
        if type(iterable[0]) not in (list, tuple):
            iterable = (iterable,) 
            
        self._thelen = len(iterable[0]) #ks init
        self._maindict = OrDi(enumerate(iterable)) #ks init
        self._orders = [] #ks init
        for i in range(self._thelen): #joka dimensiolle
            self._orders.append(([],[])) #lisätään listapari
            for key in sorted(self._maindict, key=lambda x: self._maindict[x][i]): #jokaiselle pisteelle (koordinaattien lukujärjestyksessä)
                self._orders[i][0].append(  key  ) #listaan 0 lisätään indeksi
                self._orders[i][1].append(  self._maindict[key][i]  ) #listaan 0 lisätään koordinaatti
        
    def findinrange(self, coord, radius, indexes=None , level = 0): #kutsuttaessa anna hakukoordinaatti ja säde, muut auttavat rekursiivisessa dimensioiden käsittelyssä
        '''
        tehdään ensin haut joilla haetaan säteen sisällä olevat pisteen yksittäisten dimensioiden mukaan ja vasta siten läytyneille pisteille tehdään todellinen etäisyyshaku.
        Eli ensin etsitään säteen mukainen kuutio ja vasta sen sisältä pallo :P (kun kolme ulottuvuutta)
        '''
        
        if not indexes:
            indexes = self._maindict.keys() #haetaan alkuperäiset indeksit
        
        ind_left  = bi.bisect_right( self._orders[level][1] , coord[level]-radius ) #puolitushaku piste - radius
        ind_right =  bi.bisect_left( self._orders[level][1] , coord[level]+radius ) #puolitushaku piste + radius
        
        if indexes is not self._maindict.keys(): #mikäli ei ensimmäinen taso
            new_indexes = list( set(indexes) & set(self._orders[level][0][ind_left:ind_right]) ) #ulosmenevistä indekseistä poistetaan ne joita ei ollut edellisessä haussa
        else:
            new_indexes = self._orders[level][0][ind_left:ind_right] #indeksit puolitushakujen välistä läytyvät
        
        if len(new_indexes)==0: #mikäli ei läytynyt mitään, palaa
            return new_indexes
        elif level<self._thelen-1: #mikäli ei olla viimeisessa dimensiossa
            return self.findinrange(coord, radius, new_indexes, level +1) #haetaan myäs seuraavasta dimensiosta, annetaan myäs täällä läytyneet
        else: #mikäli ollaan viimeisessä dimensiossa
            final = []
            co1 = numpy.array(coord)
            for i in new_indexes: #käydään läpi läydetyt pisteet
                if  (sum((co1-numpy.array(self._maindict[i]))**2))**0.5 <= radius: #etäisyystesti
                    final.append(i)
            return final

