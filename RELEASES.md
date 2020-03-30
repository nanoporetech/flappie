#  Flappie

Flappie is a prototype basecaller for the sequence platforms produced by Oxford Nanopore Technologies.

##  2.1.3 (_Pelotretis flavilatus_)
New Flip-flop model for direct RNA
- Accuracy improvement, typically +1.5% over Guppy Hac v3.5.1

##  2.1.2 (_Pleuronichthys decurrens_)
Further improvements to Flip-flop and RLE models
- Accuracy improvement, typically +1.5% over Guppy Hac v3.4
- Performance reduced

##  2.1.1 (_Solea solea_)
Improved Flip-flop and RLE models
- Accuracy improvements, typically > +1%
- Performance improvements, approximately 2x faster

##  2.1.0 (_Lepidopsetta bilineata_)
_Runnie(r)_, improved Run-Length Encoding.
- Improved RLE model
- Shift-invariant formulation; input reads nolonger need to be scaled
- Shape and scale parameters for run-distributions now decimal-encoded floats

##  2.0.0 (_Psettichthys melanostictus_)
_Runnie_, Run-Length Encoded basecaller

##  1.1.0 (_Isopsetta isolepis_)
**Previous release of 5mC had incorrect model - please upgrade**
- Simple trace viewer with support for modified bases
- Correct implementation of 5mC calling

##  1.0.0 (_Microstomus kitt_)
Public release of _Flappie_ under the terms of the Oxford Nanopore Technologies Ltd. Public Licence v1.0.

##  0.1.0 (_Parophrys vetulus_)
Initial release of _Flappie_
- Flip-flop basecalling of Fast5
- Calling of 5mC as a separate base in CpG contexts
