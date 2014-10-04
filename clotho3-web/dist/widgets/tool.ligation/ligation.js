$clotho.extensions.controller('clothoTool_ligation', function($scope, Clotho, PCR) {

	//todo - reduce reliance on these modules

	$scope.ligate_demoSets = [
		['aaaaaaaaaaA^CATG_', '^CATG_Tttggttggttgg'],
		['aaaaaaaaaaA^CATG_', 'ccaaccaaccaaA^CATG_'],
		['^CATG_Ttttttttttt', '^CATG_Tttggttggttgg'],
		['^CATG_Ttttttttttt', 'ccaaccaaccaaA^CATG_'],
		['aactgatcgaA^CATG_', 'acgactaA^CATG_'],
		['gctgctagcA^CATG_', 'gatcgatacc_GTAC^'],
		['acgttgcA^CATG_T', 'A^CATG_Tcagctgatgcgtcgac'],
		['ggttccggttcc|', '|tagtagtagtagtag']
	];

	$scope.setLigate = function(setInd) {
		$scope.fragments = $scope.ligate_demoSets[setInd];
	};

	function wrapFragments (fragments) {
		return angular.map(fragments, function (frag) {
			return {
				sequence : frag
			}
		});
	}

	$scope.$watch(function() {
		return $scope.fragments[0] + $scope.fragments[1]
	}, function (newval, oldval) {
		Clotho.run('clotho.functions.dna.ligate', [wrapFragments($scope.fragments)], {mute : true})
		.then(function (result) {
			$scope.ligated = result;
		});
	});

	//init
	$scope.setLigate(0);
});