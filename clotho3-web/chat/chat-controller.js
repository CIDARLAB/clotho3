'use strict';

Application.Chat.controller('ChatCtrl', ['$scope', 'Clotho', 'ChatSocket', '$dialog', '$keypress', function($scope, Clotho, ChatSocket, $dialog, $keypress) {

    /*** chat ***/

    $scope.messages = ChatSocket.messages;
    $scope.users = ChatSocket.users;

    $scope.chat_send = function(msg) {
        ChatSocket.chat_send(msg);
    };

    $scope.$on('$destroy', Clotho.silence($scope.$id));

    /***** modal ******/

    $scope.modal = {};

    $scope.modal.open = function () {
        $scope.modal.shouldBeOpen = true;
    };

    $scope.modal.close = function () {
        $scope.modal.closeMsg = 'I was closed at: ' + new Date();
        $scope.modal.shouldBeOpen = false;
    };

    $scope.modal.items = ['item1', 'item2'];

    $scope.modal.opts = {
        backdropFade: true,
        dialogFade:true
    };

    /******** dialog ********/

    // Inlined template for demo
    var dialog_template = '<div class="modal-header">'+
        '<h3>This is the title</h3>'+
        '</div>'+
        '<div class="modal-body">'+
        '<p>Enter a value to pass to <code>close</code> as the result: <input ng-model="result" type="text" /></p>'+
        '<p>{{"first letter should be capitalized - testing filter" | capitalize}}</p>' +
        '</div>'+
        '<div class="modal-footer">'+
        '<button ng-click="close(result)" class="btn btn-primary" >Close</button>'+
        '</div>';

    $scope.dialog = {};

    $scope.dialog.opts = {
        backdrop: true,
        keyboard: true,
        backdropClick: true,
        template:  dialog_template, // OR: templateUrl: 'path/to/view.html',
        controller: 'TestDialogController',
        dependencies : [
            "partials/DialogTestCtrl.js",
            "extensions/capitalize-filter.js"
        ]
    };

    $scope.openDialog = function(){
        var d = $dialog.dialog($scope.dialog.opts);
        d.open().then(function(result){
            if(result)
            {
                console.log('dialog closed with result: ' + result);
            }
        });
    };

    $scope.openMessageBox = function(){
        var title = 'This is a message box';
        var msg = 'This is the content of the message box';
        var btns = [{result:'cancel', label: 'Cancel'}, {result:'ok', label: 'OK', cssClass: 'btn-primary'}];

        $dialog.messageBox(title, msg, btns)
            .open()
            .then(function(result){
                console.log('dialog closed with result: ' + result);
            });
    };

    $keypress.on('keypress', {'enter' : 'foo()'}, $scope);

    $scope.foo = function() {
        console.log("blah");
    }
}]);