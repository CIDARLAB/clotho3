$clotho.extensions.controller('clothoTool_pcr', function($scope, DNA) {

	$scope.DNA = DNA;

	$scope.pcr_demoSets = [
		{
			primers : ['tatcgatcgta', 'gatcgatcgat'],
			backbone : 'ctcgatagcTACGATCGATAacgtagatcgatcagctgctagctagctactgaGATCGATCGATacgtacgtagctacgtacgac'
		},
		{
			primers : ['CGGATCGATCGATCGT', 'GATCGATCGATAC'],
			backbone : 'actgactacgatcgatACGATCGATCGATCCGatcgacgtaggacGATCGATCGATACacgatcgatcagctacgatcgatcgac'
		},
		{
			primers : ['GTAGTAGTAGTAGTA', 'TTCGATCGATCGATCGA'],
			backbone : 'acgtacgtacgTACTACTACTACTACagctagctacgatcgatcgatcgagcatctagctatcgtagctagcgactacgtacgtacgtcgatcgatcgatcgatcgatcgatcgactgatcgactgactgactagctacgatctacgtagctTCGATCGATCGATCGAtagctagctagctagctcga'
		},
		{
			primers : ['GTAGTAGTAA', 'TCGATCGAA'],
			backbone : 'gatcgatcgactgactagcatcgatcgatcgatcgatcgatcgagctacactgactactactattacTACTACTACTACTACacgtacgactacgacTCGATCGATCGATCGAactactacgatcgatcgatcgatcatc'
		}
	];

	$scope.setPCR = function(setInd) {
		var pcrSet = $scope.pcr_demoSets[setInd];
		$scope.primers = pcrSet.primers;
		$scope.backbone = pcrSet.backbone;
	};

	$scope.setPCR(0);
});