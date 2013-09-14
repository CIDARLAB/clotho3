'use strict';

Application.Extensions.controller('clothoIntro_toolsIntroCtrl', ['$scope', '$focus', '$timeout', '$dialog', 'Clotho', function($scope, $focus, $timeout, $dialog, Clotho) {

    $scope.showMeHow = function() {
        $focus.typeOut($('#seq'), $scope.randomSequence, 'sequence')
        .then(function() {
            return $focus.highlightElement($('#clothorun'))
        })
        .then(function(unhighlight) {
            return $timeout(function() {unhighlight() }, 1000);
        })
    };

    $scope.randomSequence = 'ACGTACATCGCGAT'; //todo - random

    $scope.gradeSimpleQuiz = function() {
        Clotho.gradeQuiz($scope.randomSequence, $scope.quizAnswer, 'aa7f191e810c19729de86101').then(function (r) {
            console.log(r);
            $scope.submitted = true;
            $scope.response = {};
            $scope.response.result = r;
        })
    };

}]);