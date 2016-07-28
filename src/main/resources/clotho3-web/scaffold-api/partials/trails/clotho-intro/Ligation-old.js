'use strict';

$clotho.extensions.controller('clothoIntro_LigationCtrl', function($scope, $focus, $timeout, $modal, DNA, Digest, PCR) {
    $scope.DNA = DNA;
    $scope.Digest = Digest;
    $scope.PCR = PCR;

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

    $scope.setLigate(0);

    $scope.$watch(function() {
        return $scope.fragments[0] + $scope.fragments[1]
    }, function (newval, oldval) {
        console.log(newval);
        $scope.ligated = PCR.ligate($scope.fragments);
    });


    $scope.showMeHow = function() {

        //todo - move to funciton that takes array to type out commands

        $modal.messageBox('Defining Variables', 'This time we\'ll define the variables as part of the process. First we need to define our two fragments, <code>frag1</code> and <code>frag2</code>. Then we\'ll join them in an array <code>fragments</code>, which we pass to the function <code>ligate</code>.', [{label: "OK", cssClass: "btn-primary", result: true}])
	          .result
            .then(function() {
                var str = 'var frag1 = "aaatttcccgggA^CATG_";';
                return $focus.typeOutSearch(str, true)
            })
            .then(function() {
                var str = 'var frag2 = "^CATG_Tttggttggttgg";';
                return $focus.typeOutSearch(str, true)
            })
            .then(function() {
                var str = 'var fragments = [frag1, frag2];';
                return $focus.typeOutSearch(str, true)
            })
            .then(function() {
                return $focus.typeOutSearch('PCR.ligate(fragments)', true);
            })
    };

});