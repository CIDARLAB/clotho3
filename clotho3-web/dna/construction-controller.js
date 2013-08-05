'use strict';

Application.Dna.controller('constructionCtrl', ['$scope', 'DNA', 'Digest', 'PCR', function($scope, DNA, Digest, PCR) {
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


    //PCR predict
    $scope.primers = ['tatcgatcgta', 'gatcgatcgat'];
    $scope.backbone = 'cccccccccccagctacgatcgataaaaaaaaaaattttttttttttgatcgatcgatagctaggggggggggggg';

}]);