/**
 * @name ClothoCommandHistory
 *
 * @description
 * Service which retrieves prior user actions, and tracks interaction with the server + and messages on say() in session
 *
 * Used by terminal and activity log
 *
 * //todo - move outside core and decorate API or PubSub? Hard because each function should post different things
 *
 * //todo - maybe create PubSub listener + publisher
 */
angular.module('clotho.core')
.service('ClothoCommandHistory', function(PubSub, $filter) {

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

		/*
		expects a well-formed message e.g. from ClientAPI.say(), of form:

		 {
			 "text" : "Welcome to Clotho!",
			 "from" : "server",
			 "class" : "success",
			 "channel" : "submit"
			 "timestamp" : Date.now()
		 }
		*/
		function addMessage (data) {
			data = angular.extend({
				"timestamp" : Date.now()
			}, data);
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

		return {
			entries : entries,

			addMessage : addMessage,

			getLastView : function () {
				return lastView;
			},
			setLastView : function (date) {
				lastView = date || Date.now();
				setUnreadCount(lastView);
			},

			toggleTerminal : function () {
				return PubSub.trigger('toggleTerminalActive');
			}
		}

});
