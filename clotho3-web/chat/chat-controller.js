'use strict';

Application.Chat.controller('ChatCtrl', ['$scope', 'Clotho', 'ChatSocket', function($scope, Clotho, ChatSocket) {

    //placeholder - model request
    $scope.uuid = 'inst_first';

    //chat
    $scope.messages = ChatSocket.messages;
    $scope.users = ChatSocket.users;

    /*** basic events ***/

    /*testing
    $scope.$watch('uuid', function(newVal, oldVal) {
        //note - will change here, not in directive (proto inheritance)
        if (newVal !== oldVal) {
        console.log("CHATCTRL\tuuid changed: " + newVal);
        }
    });
    */

    $scope.request_getModel = function(uuid) {
        $scope.requested = Clotho.get(uuid);
    };

    /*** chat ***/

    $scope.chat_send = function(msg) {
        ChatSocket.chat_send(msg);
    };

    $scope.$on('$destroy', Clotho.silence($scope.$id));

}]);