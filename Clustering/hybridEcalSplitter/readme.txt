

implementation of Kotera's strip split method, extended to hybrid case

Daniel Jeans, 2/11/2010


example steering snippet:

 <processor name="MyHybridSplitter" type="hybridRecoProcessor">
   <!--- input collections -->
   <parameter name="ECALcollections_cells"      type="StringVec"> ECALSiBarrel ECALSiEndcap </parameter>
   <parameter name="ECALcollections_tranStrips" type="StringVec"> ECALScTransverseBarrel ECALScTransverseEndcap </parameter>
   <parameter name="ECALcollections_longStrips" type="StringVec"> ECALScLongitudinalBarrel ECALScLongitudinalEndcap </parameter>
 </processor>

