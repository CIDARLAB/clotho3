/**
 * @name ClothoCommandHistory
 *
 * @description
 * Service which retrieves prior user actions, and tracks interaction with the server + and messages on say() in session
 *
 * Used by terminal and activity log
 *
 * //todo - CommandBar should use this
 * //todo - handle publishing on PubSub (clotho.trigger)
 * //todo - handle subscribing so can easily post to it (clotho.listen)
 */
angular.module('clotho.foundation')
.service('ClothoCommandHistory', function(Clotho, $filter) {

		var entries = [],
			lastView = Date.now();

			entries.unread = 0;

		//given a date, returns array of messages which are more recent
		function findNewerMessages (date) {
			date = date || lastView;
			return $filter('filter')(entries, function (item) {
				return item.timestamp > date;
			});
		}

		function setUnreadCount (date) {
			date = date || lastView;
			entries.unread = (findNewerMessages(date)).length;
		}

		//expects a well-formed message e.g. from Clotho.say()
		function addMessage (data) {
			entries.unshift(data);
			setUnreadCount();
		}

		//init
		addMessage({
			"text" : "Welcome to Clotho!",
			"from" : "server",
			"class" : "success",
			"timestamp" : Date.now()
		});

		//listeners

		//todo - change event name
		Clotho.listen("activityLog", addMessage, 'ClothoCommandHistory');

		return {
			entries : entries,

			getLastView : function () {
				return lastView;
			},
			setLastView : function (date) {
				lastView = date || Date.now();
				setUnreadCount(lastView);
			},

			toggleTerminal : function () {
				return Clotho.trigger('toggleTerminalActive');
			}
		}

});
