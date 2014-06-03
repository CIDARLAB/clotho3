$clotho.extensions.controller('clothoTool_pcr', function($scope, Clotho, DNA) {

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

})
.directive('pcrToolPredict', function (Clotho) {
	return {
		restrict: 'A',
		scope: {
			backbone: '=',
			primers: '='
		},
		link: function pcrPredictLink(scope, element, attrs) {

			function wrapSeq (string) {
				return {sequence : string};
			}

			function process () {
				Clotho.run('clotho.functions.dna.pcr', [wrapSeq(scope.backbone), wrapSeq(scope.primers[0]), wrapSeq(scope.primers[1])], {mute : true})
				.then(function (result) {
					scope.predicted = result;
					element.text(result.sequence);
				});
			}

			scope.$watch(function () {
				return scope.primers[0] + scope.primers[1] + scope.backbone;
			}, function () {
				process();
			});
		}
	};
});