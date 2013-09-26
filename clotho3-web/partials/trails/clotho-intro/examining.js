'use strict';

Application.Extensions.controller('clothoIntro_ExaminingCtrl', ['$scope', '$focus', '$timeout', '$dialog', '$http', 'DNA', 'Digest', 'PCR', 'Construction', 'Clotho', function($scope, $focus, $timeout, $dialog, $http, DNA, Digest, PCR, Construction, Clotho) {

    $http.get('/models/construction_gfp.json').then(function(data) {
        $scope.constructionFile = data.data
    });

    $scope.reactions = Construction.reactions;

    $scope.hideOpts = {
        dictionary : true,
        addStep: true,
        uptostep : 0
    };

}]);