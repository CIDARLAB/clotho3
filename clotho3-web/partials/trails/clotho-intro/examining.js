'use strict';

Application.Extensions.controller('clothoIntro_ExaminingCtrl', ['$scope', '$focus', '$timeout', '$dialog', '$http', 'DNA', 'Digest', 'PCR', 'Construction', 'Clotho', function($scope, $focus, $timeout, $dialog, $http, DNA, Digest, PCR, Construction, Clotho) {

    $scope.constructionFile = $http.get('/models/construction_gfp.json').then(function(data) { return data.data });

    $scope.$watch('constructionFile', function (newval) {
        if (!newval) return;
        $scope.constructionFileProduct = Construction.process(newval)
    });

}]);