'use strict';

Application.Extensions.controller('clothoIntro_PCRCtrl', ['$scope', '$focus', '$timeout', '$dialog', 'DNA', 'Digest', 'PCR', function($scope, $focus, $timeout, $dialog, DNA, Digest, PCR) {

    $scope.pcr_demoSets = [
        {
            primers : ['tatcgatcgta', 'gatcgatcgat'],
            backbone : 'cccccccccccagctacgatcgataaaaaaaaaaattttttttttttgatcgatcgatagctaggggggggggggg'
        },
        {
            primers : ['CGGATCGATCGATCGT', 'GATCGATCGATAC'],
            backbone : 'tttttttttttttttttACGATCGATCGATCCGggggggggggggggggggcccccccccccccccGATCGATCGATACaaaaaaaaaaaaaa'
        },
        {
            primers : ['GTAGTAGTAGTAGTA', 'TCGATCGATCGATCGAA'],
            backbone : 'ttttttttttttttttttttttttttttttttttttttTACTACTACTACTACgggggggggggggggggggggggggggggggggggggggggggggggggggggggggccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccTCGATCGATCGATCGAaaaaaaaaaaaaaaaaaaaaa'
        },
        {
            primers : ['GTAGTAGTAA', 'TCGATCGAAA'],
            backbone : 'ttttttttttttttttttttttttttttttttttttttTACTACTACTACTACgggggggggggggggggggggggggggggggccccccccccccccccccccccccccccccTCGATCGATCGATCGAaaaaaaaaaaaaaaaaaaaaa'
        }
    ];

    $scope.setPCR = function(setInd) {
        var pcrSet = $scope.pcr_demoSets[setInd];
        $scope.primers = pcrSet.primers;
        $scope.backbone = pcrSet.backbone;
    };

    $scope.setPCR(2);

    $scope.showMeHow = function() {

        //todo - dialog to DIY

        var str = 'clotho.run("PCR", [<backbone>, <primers>]';
        $focus.typeOutSearch(str)
        .then(function() {
            $('#searchBarInput').focus()
        });
    };

}]);