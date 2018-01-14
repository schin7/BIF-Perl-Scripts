use strict;
use warnings;

#Initialize variables
#Hash of known amino acids molecular weights
my %molWeights = (
	'89.0929' => "Ala", 
	'121.1579' => "Cys", 
	'133.1024' => "Asp", 
	'147.1289' => "Glu",
	'165.1887' => "Phe",
	'75.0664' => "His",
	'155.1542' => "Ile",
	'146.1870' => "Lys",
	'131.1724' => "Leu",
	'149.2109' => "Met",
	'132.1176' => "Asn",
	'115.1301' => "Pro",
	'146.1441' => "Gln",
	'174.2004' => "Ser",
	'119.1188' => "Thr",
	'117.1459' => "Val",
	'204.2247' => "Trp",
	'181.1881' => "Tyr"
); 
#Array of amino acid weights
my @AAarray = (
	'89.0929', 
	'133.1024', 
	'147.1289',
	'165.1887',
	'75.0664',
	'155.1542',
	'146.1870',
	'131.1724',
	'132.1176',
	'115.1301',
	'146.1441',
	'174.2004',
	'119.1188',
	'117.1459',
	'204.2247',
	'181.1881'
); 
#Hash of atom valences 
my %Valence = (
	"H" => 1, 
	"O" => (-2), 
	"N" => (-3), 
	"C" => (-4)
); 
#Hash of atoms 
my %AtomOutput = (
	1 => "H", 
	2 => "O", 
	3 => "N", 
	4 => "C"
); 
#Hash of atom molecular weights
my %AtomWeight = (
	"H" => 1.0079, 
	"O" => 15.9994, 
	"N" => 14.0067, 
	"C" => 12.0107
); 

my $BV = 0; 		 #Bond Value of molecule
my $BVVal = 1;
my $CarbCounter = 0; #counts carbons for occurring double bonds
my $MW = 0;          #Molecular Weight of Molecule created
my $AtomInput;       #Random number between 1-4
my $CurrentAtom; 
my $Counter = 0;     #The number of cycles
my $CountEvery10 = 0; 
my $Molecule;        #string of the atoms in the molecule
my %ClStore = ();
my %ClStoreCurrent = ();
my $ClWeight = 9999; #Stores the weight of the closest amino acid
my $ClAA;            #Closest amino acid
my $ClCounter;       #Stores the iteration where the closest molecule appears
my $i;             
my $StableMolCounter;
my $K1;
my $K2;
my $TrackedW = 9999;
my $TrackedAA;
my $TrackedWD = 9999;
my $TrackedAAD;


until (exists($molWeights{$MW})){
	$Counter++;
	$CountEvery10++;

	while($MW < '204.2247' and $BVVal){ 		#while molecule weight is greater than trp and is stable
		$AtomInput = int( 1 + rand(4)); 		#generate random number between 1-4
		$CurrentAtom = $AtomOutput{$AtomInput}; #retrieves the atom associated
		$Molecule = $Molecule . $CurrentAtom;   #stores the atom
		$MW = $MW + $AtomWeight{$CurrentAtom};  #retrieves the Molecular Weight

		if($CurrentAtom eq "C"){ #Carbon counter
			$CarbCounter++;
		}else{
			$CarbCounter = 0;
		}
		
		if (length($Molecule) == 1){ #Initializes BV 
			$BV = $Valence{$CurrentAtom};
		}elsif($BV < 0 and $CarbCounter < 3){
			if ($CurrentAtom eq "O" or $CurrentAtom eq "N" or $CurrentAtom eq "C"){
				$BV = $BV + $Valence{$CurrentAtom} + 2;
			}else{
				$BV = $BV + 1; #Hydrogen
			}
		}elsif($BV > 0 and $CarbCounter < 3){ #BV is positive
			$BV = $BV + $Valence{$CurrentAtom};
		}

		if ($CarbCounter == 3){ #Conditions for 3 carbons
			if ($BV > 0){
				$BV = $BV - 2;
			}#other condition is >0 and there will be no change
			$CarbCounter = 0;
		}

#adds 3 hydrogens when BV < -4
		if ($BV < -4){ 
			$MW = $MW + 3.0237;
			$Molecule = $Molecule . "HHH"; 
			$BV = $BV + 3;		
		}
#Exit condition
		if($BV == 0){ 
			$BVVal = 0;
			$StableMolCounter++;
		}
	}

#stores the closest molecule
	for($i=0; $i < scalar(@AAarray); $i++){ 
		if($ClWeight> abs($AAarray[$i] - $MW)){
			$ClWeight = abs($AAarray[$i] - $MW);
			$ClAA = $molWeights{@AAarray[$i]};
			$ClCounter = $Counter;
		}
	}

#if it is not in the hash, then it will store
	if(not exists($ClStore{$ClAA})){ 
		$ClStore{$ClAA} = $ClWeight;
	}else{
		if($ClWeight < $ClStore{$ClAA}){
			$ClStore{$ClAA} = $ClWeight;	
		}
	}

	if(not exists($ClStoreCurrent{$ClAA})){ 
		$ClStoreCurrent{$ClAA} = $ClWeight;
	}else{
		if($ClWeight < $ClStoreCurrent{$ClAA}){
			$ClStoreCurrent{$ClAA} = $ClWeight;
		}
	} 
	
	if (not exists($molWeights{$MW})){
		print "Molecule: $Counter has been produced!\n";
		print "Closest amino acid (by weight): $ClAA\n";
		
		foreach $K1(keys(%ClStoreCurrent)){
			if($TrackedWD > $ClStoreCurrent{$K1}){
				$TrackedWD = $ClStoreCurrent{$K1};					
				$TrackedAAD = $K1;
			}
		}

		print "Current closest matching weight: $TrackedAAD Difference: $TrackedWD \n";

		if($CountEvery10 == 10){
			foreach $K2(keys(%ClStore)){
				if($TrackedW > $ClStore{$K2}){
					$TrackedW = $ClStore{$K2};
					$TrackedAA = $K2;
				}
			}
			print "Closest match in 10 cycles: $TrackedAA Difference: $TrackedW \n";

			$CountEvery10 = 0;
			%ClStore = ();
			$TrackedW = 9999;
			$TrackedWD = 9999;
		}
		
		#Reset program
		$Molecule = "";
		$MW = 0;
		$BVVal = 1;
		$BV = 0;
		$ClWeight = 9999;

	}else{
		print "The simulation has produced: ", $molWeights{$MW}, "\n";
		print "Stable molecules created: $StableMolCounter \n";
		print "End of Program. ";
		exit
	}
}

