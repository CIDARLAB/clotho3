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
    $scope.digestFrags = Digest.digest($scope.digestSeq, $scope.digestEnz, false);

    $scope.$watch('digestEnz', function() {
        $scope.digestFrags = Digest.digest($scope.digestSeq, $scope.digestEnz, false);
    });


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
    
    
    
    console.log(PCR.findAnnealSimple($scope.backbone, $scope.primers[0]));


    //ligation

    //first 4 should all be the same
    $scope.ligate_demoSets = [
        ['aaaaaaaaaaA^CATG_', '^CATG_Tttggttggttgg'],
        ['aaaaaaaaaaA^CATG_', 'ccaaccaaccaaA^CATG_'],
        ['^CATG_Ttttttttttt', '^CATG_Tttggttggttgg'],
        ['^CATG_Ttttttttttt', 'ccaaccaaccaaA^CATG_'],
        ['aaaaaaaaaaA^CATG_', 'ggggggA^CATG_'],
        ['aaaaaaaaaaA^CATG_', 'gtcatcgatcagt_GTAC^'],
        ['aaaaaaaaaaA^CATG_T', 'A^CATG_Tacgatagcattaagcgt'],
        ['_gggg^tttttttttttt', 'aaaaa^cccc_'],
        ['aaaaaaaaaaaaaa|', '|ttggttggttgg']
    ];

    $scope.setLigate = function(setInd) {
        $scope.fragments = $scope.ligate_demoSets[setInd];
    };

    $scope.setLigate(0);


    $scope.$watch('fragments', function () {
        $scope.ligated = PCR.ligate($scope.fragments);
    }, true);


    $scope.clothoFunctions = [];
    Clotho.query({schema: "Function"}).then(function(result) {
        console.log(result);
        $scope.clothoFunctions = result;
    });


}]);