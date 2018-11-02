
## ************************************************************************
## 
## 
## 
## (c) Xiaobei Zhao
##
## Fri Nov 02 00:59:54 EDT 2018 -0400 (Week 44)
## Fri Mar 17 13:41:53 EDT 2017 -0400 (Week 11)
## 
## 
## Reference: 
##
## [Note]
##   This is to convert a data table in a text file to VCF format.
##   The columns are rearranged according to VCF's specifications;
##   the content is kept intact.
## [README]
## To run an example:
##   python convert.py
##
## 
## 
## ************************************************************************


import os
import StringIO


def is_numeric(x,complex_allow=False):
  '''
  @examples
  is_numeric(None)
  '''
  ret = True
  try:
    if complex_allow:
      _x=complex(x)
    else:
      _x=float(x)
    ##
    ret=(_x == _x)
  except Exception: #ValueError, TypeError
    ret=False
  ##
  return(ret)

  
def get_colnames(
  inFpath,
  header=True,
  sep="\t",
  skip=None,
  close=True
  ):
  
  if skip is None:
    skip=0

  ## i/o
  if type(inFpath) in [type(StringIO.StringIO(""))]:
    _input=inFpath
  else:
    _input=open(inFpath)

  ##
  if header:
    for i in xrange(skip):
      next(_input)
    ##
    row=_input.readline().strip()
    _row=row.rstrip("\n").split(sep)
    ret=_row
  else:
    ret=None

  ##
  if close and not _input.closed:
    _input.close()
  ##
  return(ret)
      


def txt_to_vcf(
  inFpath,
  outFpath=None,
  header=True,

  sep_in="\t",
  sep_out="\t",
  skip=None,
  colnames_in=None,
  colnames_out=["CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO"], #FORMAT (TBA)
  colnames_mapping={"CHROM":1,"POS":2,"ID":None,"REF":3,"ALT":4}, #or 1-baesd index

  info_name="INFO",
  info_sep=";",
  info_equ="="
  ):

  if sep_in is None:
    sys.exit("`sep_in` must be specified.")

  if skip is None:
    skip=0

  ##
  if colnames_in is None:
    _colnames=get_colnames(inFpath=inFpath,header=header,sep=sep_in,skip=skip,close=False)
    if header:
      colnames_in=_colnames
    else:
      colnames_in=["V"+str(i+1) for i in xrange(len(_colnames))]
  
  ## check names
  # # print(colnames_in)
  # # print(colnames_mapping)
  _colnames_mapping={k:colnames_in.index(v)+1  if (not is_numeric(v) and not v is None) else colnames_mapping[k] for (k,v) in colnames_mapping.iteritems()} # using 1-based index instead of col_name
  _values=_colnames_mapping.values()
  _whether=[not is_numeric(e) or not e is None for e in _values]
  if not all(_whether):
    _msg=",".join([a for (a,b) in zip(_values,_whether) if not b])
    sys.exit("Unrecognized column names: "+_msg)
  _colnames_mapping={k:None  if (v is None or v > len(colnames_in)) else _colnames_mapping[k] for (k,v) in _colnames_mapping.iteritems()}
  _colnames_mapping={k:v-1  if (v is not None) else _colnames_mapping[k] for (k,v) in _colnames_mapping.iteritems()}
  # # print(_colnames_mapping)
  
  if info_name in colnames_out:
    info_index=colnames_out.index(info_name)
  else:
    info_index=None
    
  ## i/o
  if type(inFpath) in [type(StringIO.StringIO(""))]:
    _input=inFpath
  else:
    _input=open(inFpath)

  ##
  if outFpath is None:
    if type(inFpath)==type(""):
      outFpath=inFpath+".vcf"
    
  if outFpath is None:
    _output=StringIO.StringIO()
  else:
    _output=open(outFpath,'w')
    
    
  ## write before_header
  ## (TBA)
  ## write header
  _output.write("#%s\n" % sep_out.join(colnames_out))

  ##
  _input.seek(0)
  ##
  for i in xrange(skip):
    next(_input)
  ##
  if header:
    next(_input)
  ##
  while True:
    row=_input.readline().strip()
    if not row:
      break
    ##
    _row=row.split(sep_in)
    tmp=[_row[_colnames_mapping[k]] if k in _colnames_mapping and _colnames_mapping[k] is not None else "." for k in colnames_out]
    ##
    if not info_index is None:
      _info=[(colnames_in[i],_row[i]) for i in xrange(len(_row)) if i not in _colnames_mapping.values()]
      _info=info_sep.join([info_equ.join((a,b)) for (a,b) in _info])
      ##
      if tmp[info_index]==".":
        tmp[info_index]=_info
      else:
        tmp[info_index]=info_sep.join([tmp[info_index],_info])
    ##
    _output.write("%s\n" % sep_out.join(tmp))
  ##
  if type(_input)==type(StringIO.StringIO("")):
    print("input:")
    print(_input.getvalue())
    print("")
  if type(_output)==type(StringIO.StringIO("")):
    print("output:")
    print(_output.getvalue())
    print("")
  ##
  if not _input.closed:
    _input.close()
  ##
  if not _output.closed:
    _output.close()
  ##
    
  return()


if __name__ == "__main__":

  input=StringIO.StringIO(
    '''
id;chr;pos;ref;alt;gene
260;chr1;38180647;T;C;EPHA10
1310;chr2;103281538;G;GT;SLC9A2
2494;chr4;106164913;C;A;TET2
2593;chr4;187513886;G;A;FAT1
3034;chr5;176637352;C;T;NSD1
3103;chr6;10404781;C;T;TFAP2A
7333;chr15;41023367;C;T;RAD51
7948;chr16;66430034;G;C;CDH5
8312;chr17;40481624;A;C;STAT3
8675;chr18;48573566;A;G;SMAD4
8712;chr18;52942953;G;A;TCF4
9468;chr20;39750756;C;G;TOP1
10140;chrX;47044518;G;A;RBM10
10154;chrX;63413026;G;C;AMER1
10211;chrX;133507309;C;G;PHF6

'''.strip())
  txt_to_vcf(inFpath=input,outFpath=None,sep_in=";",sep_out="\t",colnames_mapping={"CHROM":"chr","POS":"pos","REF":"ref","ALT":"alt"})
  input.close()

  ## [Output]
'''
input:
id;chr;pos;ref;alt;gene
260;chr1;38180647;T;C;EPHA10
1310;chr2;103281538;G;GT;SLC9A2
2494;chr4;106164913;C;A;TET2
2593;chr4;187513886;G;A;FAT1
3034;chr5;176637352;C;T;NSD1
3103;chr6;10404781;C;T;TFAP2A
7333;chr15;41023367;C;T;RAD51
7948;chr16;66430034;G;C;CDH5
8312;chr17;40481624;A;C;STAT3
8675;chr18;48573566;A;G;SMAD4
8712;chr18;52942953;G;A;TCF4
9468;chr20;39750756;C;G;TOP1
10140;chrX;47044518;G;A;RBM10
10154;chrX;63413026;G;C;AMER1
10211;chrX;133507309;C;G;PHF6

output:
#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO
chr1    38180647        .       T       C       .       .       id=260;gene=EPHA10
chr2    103281538       .       G       GT      .       .       id=1310;gene=SLC9A2
chr4    106164913       .       C       A       .       .       id=2494;gene=TET2
chr4    187513886       .       G       A       .       .       id=2593;gene=FAT1
chr5    176637352       .       C       T       .       .       id=3034;gene=NSD1
chr6    10404781        .       C       T       .       .       id=3103;gene=TFAP2A
chr15   41023367        .       C       T       .       .       id=7333;gene=RAD51
chr16   66430034        .       G       C       .       .       id=7948;gene=CDH5
chr17   40481624        .       A       C       .       .       id=8312;gene=STAT3
chr18   48573566        .       A       G       .       .       id=8675;gene=SMAD4
chr18   52942953        .       G       A       .       .       id=8712;gene=TCF4
chr20   39750756        .       C       G       .       .       id=9468;gene=TOP1
chrX    47044518        .       G       A       .       .       id=10140;gene=RBM10
chrX    63413026        .       G       C       .       .       id=10154;gene=AMER1
chrX    133507309       .       C       G       .       .       id=10211;gene=PHF6

'''

