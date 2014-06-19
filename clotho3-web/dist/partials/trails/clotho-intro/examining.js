'use strict';

$clotho.extensions.controller('clothoIntro_ExaminingCtrl', function($scope, $focus, $timeout, $http, Construction, Clotho) {

    $http.get('models/construction-old/construction_gfp.json').then(function(data) {
        $scope.constructionFile = data.data
    });

    $scope.reactions = Construction.reactions;

    $scope.hideOpts = {
        dictionary : true,
        addStep: true,
        uptostep : 0,
        processto: 0
    };

    $scope.revealVideo = function () {
        $scope.showVideo = true;
    };

    $scope.hideVideo = function() {
        $scope.showVideo = false;
    };

});