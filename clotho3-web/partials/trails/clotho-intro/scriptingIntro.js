'use strict';

Application.Extensions.controller('clothoIntro_scriptingIntroCtrl', ['$scope', '$focus', '$timeout', '$dialog', 'Clotho', function($scope, $focus, $timeout, $dialog, Clotho) {

    $scope.randomSequence = 'GATCGATACGATCG'; //todo - random

    $scope.gradeSimpleQuiz = function() {
        Clotho.gradeQuiz($scope.randomSequence, $scope.quizAnswer, 'aa7f191e810c19729de86101').then(function (r) {
            console.log(r);
            $scope.submitted = true;
            $scope.response = {};
            $scope.response.result = r;
        })
    };

    $scope.showMeHow = function() {
        var searchInput = $('#searchBarInput'),
            searchSubmit = $('#searchBarSubmit'),
            oldinputZ = searchInput.css("z-index"),
            oldsubmitZ = searchSubmit.css("z-index"),
            maxZ = $focus.maxZ;


        $dialog.messageBox('Entering Commands', 'The search bar executes commands in JavaScript. You can run Clotho API commands, and functions within Clotho. For example, to reverse complement a sequence, you would call the function <code>revcomp</code> and pass your sequence in the array of arguments like so: <code>clotho.run(\'revcomp\', ["' + $scope.randomSequence + '"]</code>', [{label: "OK", cssClass: "btn-primary", result: true}]).open()
        .then(function() {
            return $focus.typeOutSearch("clotho.run('revcomp', ['"+$scope.randomSequence+"'])")
        })
        .then(function() {
            return $dialog.messageBox('Submit', 'You can submit your command by either hitting the enter key, or pressing the submit button', [{label: "OK", cssClass: "btn-primary", result: true}]).open()
        });
    };


    $scope.showHelpTips = function() {
        $focus.elementPopover('#searchBarInput', "whats up");
    };


}]);