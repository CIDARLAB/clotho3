$clotho.extensions.controller('clothoIntro_reverseComplementCtrl', function($scope, $focus, $timeout, $modal, Clotho) {

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

        $modal.open(dialogOpts)
	        .opened
	        .then(function() {
		        $focus.highlightElement($('.quiz'))
	        });

});