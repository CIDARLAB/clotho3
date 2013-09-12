'use strict';

Application.Extensions.controller('clothoIntro_LigationCtrl', ['$scope', '$focus', '$timeout', '$dialog', 'DNA', 'Digest', 'PCR', function($scope, $focus, $timeout, $dialog, DNA, Digest, PCR) {
    $scope.DNA = DNA;
    $scope.Digest = Digest;
    $scope.PCR = PCR;

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

}]);