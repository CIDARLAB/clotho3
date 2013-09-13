'use strict';

Application.Extensions.controller('clothoIntro_RestrictionDigestsCtrl', ['$scope', '$focus', '$timeout', '$dialog', 'Digest', function($scope, $focus, $timeout, $dialog, Digest) {
    $scope.Digest = Digest;

    $scope.demo = {};
    $scope.demo.digestEnz = Digest.enzymes.EcoRI;
    $scope.demo.digestSeq = 'acaacgtctcacggatccagtcggaattctacatgcatcgatcgacggatccagatcgactagc';


    $scope.showMeHow = function () {
        $focus.typeOutSearch('clotho.run("digest", [<sequence>, <enzyme>])');
    }
}]);