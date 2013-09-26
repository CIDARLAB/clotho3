'use strict';

Application.Extensions.controller('clothoIntro_toolsIntroCtrl', ['$scope', '$focus', '$timeout', '$dialog', 'Clotho', function($scope, $focus, $timeout, $dialog, Clotho) {

    var dialogOpts = {
        backdrop: false,
        keyboard: false,
        dialogFade : true,
        templateUrl: '/interface/templates/dialogMessagebox.html',
        controller: 'MessageBoxController',
        resolve:
        {model: function() {
            return {
                title: 'Tools in Clotho',
                message: 'Clotho can automate several common tasks, such as reverse complementing a sequence. Clotho provides tools - like this one - to streamline the process.',
                buttons: [{label: "OK", cssClass: "btn-primary", result: true}]
            };
        }}
    };

    var dialog2Opts = {
        backdrop: false,
        keyboard: false,
        dialogFade : true,
        templateUrl: '/interface/templates/dialogMessagebox.html',
        controller: 'MessageBoxController',
        resolve:
        {model: function() {
            return {
                title: 'Show Me How',
                message: 'Buttons like this will illustrate the process for you.',
                buttons: [{label: "OK", cssClass: "btn-primary", result: true}]
            };
        }}
    };

    //allow template to render
    $timeout(angular.noop,200).
    then(function() {
        $focus.addBackdrop();
        return $focus.bringToFront($('#demoTool'))
    })
    .then(function(oldZ) {
        return $timeout(function() {return oldZ}, 500 );
    })
    .then(function(oldZ) {
        return $dialog.dialog(dialogOpts).open().then(function() {return oldZ});
    })
    .then(function(oldZ) {
        return $timeout(function() {
            $focus.setZ(oldZ, $('#demoTool'));
            return $focus.bringToFront($('#showMeHow'));
        }, 200);
    })
    .then(function(oldZ) {
        return $dialog.dialog(dialog2Opts).open().then(function() {return oldZ});
    })
    .then(function(oldZ) {
        $focus.setZ(oldZ, $('#showMeHow'));
        return $focus.removeBackdrop();
    });



    $scope.quiz = {
        "type" : "fillin",
        "hint" : "Copy-paste the sequence into the tool below",
        "question" : "Find the reverse complement of the sequence <code>{{ questionValue }}</code> using the tool.",
        "questionValue" : "ACGTACATCGCGAT",
        "answerGenerator" : "aa7f191e810c19729de86101",
        "retry" : {
            "questionValue" : "clotho.run('randomSequence', ['13'])"
        }
    };

    $scope.showMeHow = function() {
        $focus.typeOut($('#seq'), $scope.quiz.questionValue, 'sequence')
        .then(function() {
            return $focus.highlightElement($('#clothorun'))
        })
        .then(function(unhighlight) {
            return $timeout(function() {unhighlight() }, 1000);
        })
    };
}]);