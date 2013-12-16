$clotho.extensions.controller('clothoIntro_reverseComplementCtrl', function($scope, $focus, $timeout, $dialog, Clotho) {

    var dialogOpts = {
        backdrop: false,
        keyboard: false,
        dialogFade : true,
        templateUrl: 'views/_interface/ui-custom/dialogMessagebox.html',
        controller: 'MessageBoxController',
        resolve:
        {model: function() {
            return {
                title: 'Quizzes',
                message: 'Boxes like the one below are quizzes. This Trail will ask you to answer questions about DNA manipulation. Watch the videos on each page if you need to brush up on your knowledge (click OK and peek behind this window), but make sure you understand each quiz before moving forward.',
                buttons: [{label: "OK", cssClass: "btn-primary", result: true}]
            };
        }}
    };

    //hacks to hide dialog
    var dialog;
    var oldNext = $scope.next;
    $scope.$parent.next = function() {
        console.log('called it');
        console.log(dialog);
        dialog.close();
        oldNext();
    };

    //allow tide for quiz to render
    $timeout(angular.noop,200)
    .then(function() {
         $focus.addBackdrop();
        return $focus.highlightElement($('.quiz'))
    })
    .then(function() {
        return $timeout(function() {  }, 500 );
    })
    .then(function() {
        dialog = $dialog.dialog(dialogOpts);
        return dialog.open()
    })
    .then(function() {
        return $focus.removeBackdrop();
    })

});