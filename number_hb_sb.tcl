set outfile1 [open nmb_hb_first.dat a]
set outfile2 [open nmb_hb_second.dat a]
set outfile3 [open nmb_hb_third.dat a]
set outfile4 [open nmb_hb_fourth.dat a]
set outfile5 [open nmb_hb_all.dat a]
set outfile6 [open nmb_wm_r.dat a]

set out_sb1 [open nmb_sb_1.dat a]
set out_sb2 [open nmb_sb_2.dat a]

set helix1hb [open helix1.dat a]
set helix2hb [open helix2.dat a]

set distfile25_44 [open distance_25_44.dat a]
set distfile56_75 [open distance_25_44.dat a]

set nf [molinfo top get numframes]

set sel11 [atomselect top "backbone and not sidechain and resid 14 to 25"]
set sel12 [atomselect top "backbone and not sidechain and resid 3 to 13"]
set sel21 [atomselect top "backbone and not sidechain and resid 3 to 12"]
set sel22 [atomselect top "backbone and not sidechain and resid 46 to 55"]

set sel31 [atomselect top "backbone and not sidechain and resid 46 to 55"]
set sel32 [atomselect top "backbone and not sidechain and resid 86 to 94"]

set sel41 [atomselect top "backbone and not sidechain and resid 76 to 85"]
set sel42 [atomselect top "backbone and not sidechain and resid 86 to 95"]

set sel51 [atomselect top "water and name OH2 and (within 4.5 of (resid 10 29 33 36 40 50 52 63 67 71 and not water))"]

set helix1 [atomselect top "resid 25 to 44 and not sidechain"]
set helix2 [atomselect top "resid 56 to 75 and not sidechain"]

set sel_25 [atomselect top "resid 25 and name CA"]
set sel_44 [atomselect top "resid 44 and name CA"]

set sel_56 [atomselect top "resid 56 and name CA"]
set sel_75 [atomselect top "resid 75 and name CA"]


#set w_molecules [$sel51 num]
for {set i 1} {$i<=$nf} {incr i} {
         
         puts "frame $i"
         
         $sel11 frame $i
         $sel12 frame $i
         
         $sel21 frame $i
         $sel22 frame $i
         
         $sel31 frame $i
         $sel32 frame $i
         
         $sel41 frame $i
         $sel42 frame $i
         
         $helix1 frame $i
         $helix2 frame $i
         
         #$sel51 frame $i
         
         $sel_25 frame $i 
         $sel_44 frame $i 
         
         $sel_56 frame $i 
         $sel_75 frame $i
         
                                                                                                                                                                  
                                                                                                                                                                                                                                          
         set hbcount1 [llength [lindex [measure hbonds 3.5 45 $sel11 $sel12] 0]]
         set hbcount2 [llength [lindex [measure hbonds 3.5 45 $sel21 $sel22] 0]]
         set hbcount3 [llength [lindex [measure hbonds 3.5 45 $sel31 $sel32] 0]]
         set hbcount4 [llength [lindex [measure hbonds 3.5 45 $sel41 $sel42] 0]]
         set hbcount5 [[atomselect top "water and name OH2 and (within 4.5 of (resid 10 29 33 36 40 50 52 63 67 71 and not water))" frame $i] num] 
         #set hbcount5 [$sel51 frame $i num]
                                                                                                                                                                                    
         set sb_count1 [[atomselect top "protein and (resname ASP GLU and (name O or name OD1 or name OD2 or name OE1 or name OE2)) and within 3.2 of (protein and resname ARG HIS LYS HSP and (name N or name ND1 or name ND2 or name NE1 or name NE2))" frame $i] num]
                                                                                                                                                                                                      
         set sb_count2 [[atomselect top "protein and (resname ARG HIS LYS HSP and (name N or name ND1 or name ND2 or name NE1 or name NE2)) and within 3.2 of (protein and resname ASP GLU and (name O or name OD1 or name OD2 or name OE1 or name OE2))" frame $i] num]                   
         
         set helix1_hb_count [llength [lindex [measure hbonds 3.5 60 $helix1] 0]]
         
         set helix2_hb_count [llength [lindex [measure hbonds 3.5 60 $helix2] 0]]
         
         set x1 [$sel_25 get {x}] 
         set y1 [$sel_25 get {y}]         
         set z1 [$sel_25 get {z}]
         
         set x2 [$sel_44 get {x}]
         set y2 [$sel_44 get {y}]
         set z2 [$sel_44 get {z}]

         set x3 [$sel_56 get {x}] 
         set y3 [$sel_56 get {y}]         
         set z3 [$sel_56 get {z}]
         
         set x4 [$sel_75 get {x}]
         set y4 [$sel_75 get {y}]
         set z4 [$sel_75 get {z}]                                                                                                                                                                                                                        
                                                                                                                                                                           
         set broj [expr "$hbcount1 + $hbcount2 + $hbcount3 + $hbcount4"]  
         
         set rastojanje1 [expr "sqrt(($x1-$x2)*($x1-$x2) + ($y1-$y2)*($y1-$y2) + ($z1-$z2)*($z1-$z2) )"]     
         puts $distfile25_44 "$rastojanje1"
         
         set rastojanje2 [expr "sqrt(($x3-$x4)*($x3-$x4) + ($y3-$y4)*($y3-$y4) + ($z3-$z4)*($z3-$z4) )"]
         puts $distfile56_75 "$rastojanje2"
                                                                                                                                                                                                      
                  
         puts $outfile1 "$hbcount1"
         puts $outfile2 "$hbcount2"
         puts $outfile3 "$hbcount3"
         puts $outfile4 "$hbcount4"
         
         puts $outfile5 "$broj"
                  
         puts $outfile6 "$hbcount5"
         
         puts $out_sb1  "$sb_count1"
         puts $out_sb2  "$sb_count2"
         
         puts $helix1hb $helix1_hb_count 
         puts $helix2hb $helix2_hb_count 
         
}
         
close $outfile1
close $outfile2
close $outfile3
close $outfile4
close $outfile5
close $outfile6
close $out_sb1
close $out_sb2
close $helix1hb
close $helix2hb

close $distfile25_44
close $distfile56_75


unset sel11;
unset sel12;
#$sel11 delete;
#$sel12 delete;

unset sel21;
unset sel22;
#$sel21 delete;
#$sel22 delete;

unset sel31;
unset sel32;
#$sel31 delete;
#$sel32 delete;

unset sel41;
unset sel42;
#$sel41 delete;
#$sel42 delete;

unset sel51;
#$sel51 delete;

$sb_count1 delete
$sb_count2 delete


$hbcount1 delete
$hbcount2 delete
$hbcount3 delete
$hbcount4 delete
$hbcount5 delete