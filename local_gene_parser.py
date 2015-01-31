from pyparsing import Word, Keyword, opAssoc, infixNotation 
import sys

term = Word('abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789_.')

AND = Keyword("and")
OR = Keyword("or")

class Operation(object):
    def __init__(self, tokens):
        self._tokens = tokens[0]
        self.assign()

    def assign(self):
        """
        function to copy tokens to object attributes
        """

    def __repr__(self):
        return self.__class__.__name__ + ":" + repr(self.__dict__)
    __str__ = __repr__

class BinOp(Operation):
    def assign(self):
        self.op = self._tokens[1]
        self.terms = self._tokens[0::2]
        del self._tokens

# def islistofstrings(list):
#     for item in list:
#         if ~isinstance(item, str):
#             return False
#     return True

        
class AndOp(BinOp):
    def make_iter(self):
        while len(self.terms) > 2:
            new_term = AndOp([[self.terms.pop(),'and',self.terms.pop()]])
            self.terms.append(new_term)
        item1 = self.terms[0]
        item2 = self.terms[1]
        if isinstance(item1, str):
            if isinstance(item2, str):
                yield item1 + ',' + item2
            else:
                for item2_1 in item2.make_iter():
                    yield item1 + ',' + item2_1
        else:
            if isinstance(item2, str):
                for item1_1 in item1.make_iter():
                    yield item1_1 + ',' + item2
            else:
                for item1_1 in item1.make_iter():
                    for item2_1 in item2.make_iter():
                        yield item1_1 + ',' + item2_1
        
                
class OrOp(BinOp):    
    def make_iter(self):
        for e in self.terms:
            if isinstance(e, str):
                yield e
            else:
                for s in e.make_iter():
                    yield s

expr = infixNotation(term,
    [
    (AND, 2, opAssoc.LEFT, AndOp),
    (OR, 2, opAssoc.LEFT, OrOp),
    ])

def invert(gene_logic):
    try:
        logic = expr.parseString(gene_logic)
    except:
        try:
            ## Might be an incorrect encoding - try fixing it
            gene_logic = gene_logic.decode('utf-8').replace(u'\xa0', u' ')
            gene_logic = str(gene_logic)
            logic = expr.parseString(gene_logic)
        except:
            
            print("\n\nBad logic: '{}'".format(gene_logic))
            return []
    #print logic
    if isinstance(logic[0], str):
        return logic
    else:
        return logic[0].make_iter()
    
def gene_parser(gene_logic):
    """Return a list of enzymes from a logical gene combination expression."""
    enzyme_list = []
    gene_logic = gene_logic.strip()
    if len(gene_logic) > 0:
        for item in invert(gene_logic):
            
            ## Sort items in 'item'
            genes = item.split(",")
            if type(genes) != list:
                enzyme_list.append(item)
            else:
                enzyme_list.append(",".join(sorted(genes)))
            
            enzyme_list.append(item)
    enzyme_list = list(set(enzyme_list))
    return enzyme_list

if __name__ == '__main__':
    gene_logic = '((BT_1625 or BT_4136) and ((BT_4134 and BT_4135) or (BT_1630 and BT_1631)))'
    from pprint import pprint
    #print(invert(gene_logic))
    for item in invert(gene_logic):
        ## Sort items in 'item'
        genes = item.split(",")
        if type(genes) != list:
            print item
        else:
            print ",".join(sorted(genes))
    
    