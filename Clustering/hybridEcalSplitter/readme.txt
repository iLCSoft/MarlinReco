DDStripSplitter
---------------

implementation of Kotera's strip split method, for new DD4hep based models
[hybrid case ( mix of scint and si ) not yet implemented]

Daniel Jeans, Nov 2018.

example steering snippet:

  <processor name="myDDStripSplitter" type="DDStripSplitter">
    <parameter name="ECALcollections_tranStrips"> ECalBarrelScOddCollectionRec ECalEndcapScOddCollectionRec </parameter>
    <parameter name="ECALcollections_longStrips"> ECalBarrelScEvenCollectionRec ECalEndcapScEvenCollectionRec </parameter>
    <parameter name="cellSize0"> 5. </parameter>
    <parameter name="cellSize1"> 45. </parameter>
  </processor>




hybridRecoProcessor
--------------------

implementation of Kotera's strip split method, extended to hybrid case

Daniel Jeans, 2/11/2010


example steering snippet:

 <processor name="MyHybridSplitter" type="hybridRecoProcessor">
   <!--- input collections -->
   <parameter name="ECALcollections_cells"      type="StringVec"> ECALSiBarrel ECALSiEndcap </parameter>
   <parameter name="ECALcollections_tranStrips" type="StringVec"> ECALScTransverseBarrel ECALScTransverseEndcap </parameter>
   <parameter name="ECALcollections_longStrips" type="StringVec"> ECALScLongitudinalBarrel ECALScLongitudinalEndcap </parameter>
 </processor>

