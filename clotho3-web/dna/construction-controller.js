'use strict';

Application.Dna.controller('constructionCtrl', ['$scope', 'Clotho', 'DNA', 'Digest', 'PCR', function($scope, Clotho, DNA, Digest, PCR) {
    $scope.DNA = DNA;
    $scope.PCR = PCR;
    $scope.Digest = Digest;

    //todo - if possible, assign these directly from DOM using a directive...
    $scope.$watch('sequence', function(newval, oldval) {
        $scope.rna = DNA.transcribe(newval);
    });

    $scope.$watch('rna', function (newval, oldval) {
        $scope.protein = DNA.translate(newval);
    });

    $scope.sequence = 'aaatttgggcccatgcta';


    //Digest
    $scope.digestSeq = 'acaacgtctcacggatccagtcggaattctacatgcatcgatcgacggatccagatcgactagc';
    $scope.digestEnz = Digest.enzymes.EcoRI;


    //PCR predict


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
            backbone : 'ttttttttttttttttttttttttttttttttttttttTACTACTACTACTACgggggggggggggggggggggggggggggggccccccccccccccccccccccccccccccTCGATCGATCGATCGAaaaaaaaaaaaaaaaaaaaaa'
        },
        {
            primers : ['GTAGTAGTAA', 'TCGATCGAAA'],
            backbone : 'ttttttttttttttttttttttttttttttttttttttTACTACTACTACTACgggggggggggggggggggggggggggggggccccccccccccccccccccccccccccccTCGATCGATCGATCGAaaaaaaaaaaaaaaaaaaaaa'
        }
    ];

    $scope.setPCR = function(pcrSetInd) {
        var pcrSet = $scope.pcr_demoSets[pcrSetInd];
        $scope.primers = pcrSet.primers;
        $scope.backbone = pcrSet.backbone;
    };

    $scope.setPCR(2);


    //ligation
    $scope.fragments = [
        'aaaaaaaaaaaa^cccc_',
        '_gggg^tttttttttttt'
    ];

    $scope.$watch(function() {
        return $scope.fragments[0] + $scope.fragments[1]
    }, function (newval, oldval) {
        console.log(newval);
        $scope.ligated = PCR.ligate($scope.fragments);
    });


    $scope.clothoFunctions = [];
    Clotho.query({schema: "Function"}).then(function(result) {
        console.log(result);
        $scope.clothoFunctions = result;
    });


}]);