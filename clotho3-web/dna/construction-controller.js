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
    $scope.primers = ['tatcgatcgta', 'gatcgatcgat'];
    $scope.backbone = 'cccccccccccagctacgatcgataaaaaaaaaaattttttttttttgatcgatcgatagctaggggggggggggg';

    $scope.clothoFunctions = [];
    Clotho.query({schema: "Function"}).then(function(result) {
        console.log(result);
        $scope.clothoFunctions = result;
    });


}]);