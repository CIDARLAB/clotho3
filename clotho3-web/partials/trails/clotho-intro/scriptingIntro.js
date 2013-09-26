'use strict';

Application.Extensions.controller('clothoIntro_scriptingIntroCtrl', ['$scope', '$focus', '$timeout', '$dialog', 'Clotho', function($scope, $focus, $timeout, $dialog, Clotho) {

    $scope.quiz = {
        "type" : "fillin",
        "question" : "Find the reverse complement of the sequence <code>{{ questionValue }}</code> using the Command Bar.",
        "questionValue" : "GATCGATACGATCG",
        "answerGenerator" : "aa7f191e810c19729de86101",
        "retry" : {
            "questionValue" : "clotho.run('randomSequence', ['16'])"
        }
    };

    $scope.showMeHow = function() {

        $dialog.messageBox('Entering Commands', 'To reverse complement a sequence, you would call the function <code>DNA.revcomp</code> and pass your sequence in the array of arguments like so: <code>DNA.revcomp("' + $scope.quiz.questionValue + '")</code>', [{label: "OK", cssClass: "btn-primary", result: true}]).open()
        .then(function() {
            return $focus.typeOutSearch("DNA.revcomp('"+$scope.quiz.questionValue+"')")
        })
    };

}]);