//note - written using lodash internally

angular.module('clotho.core').service('PubSub',
	function ($window, $rootScope, $filter, Debug) {

		//see if already exists, don't re-instantiate
		return ($window.$clotho.$pubsub) ? $window.$clotho.$pubsub : $window.$clotho.$pubsub = generatePubSubObject();

		function generatePubSubObject() {
			/*****
			 CONFIG
			 *******/

			/*format:

			{
				topic : [
					{
						callback : <function>
						once : <boolean>
						ref : <ref>,
						priority : <number>
					}
				],
				...
			}
			*/
			var map = {};

			//split events passed in by a space
			var eventSplitter = /\s+/;

			var Debugger = new Debug('PubSub', '#55cccc');

			/*********
			 Internal Helpers
			 **********/

			//convert space-separated topic list into array
			function splitTopics (topic) {
				if (!topic) {
					return [];
				}
				else if (topic == 'all') {
					return Object.keys(map);
				}
				else {
					return topic.split(eventSplitter);
				}
			}

			function checkTopicExists(topic) {
				return map.hasOwnProperty(topic);
			}

			function checkTopicHasSubs(topic) {
				return checkTopicExists(topic) && map[topic].length > 0;
			}

			//create subscriber object
			function createSubscriber (callback, ref, priority, once) {
				if (!angular.isFunction(callback)) {
					return null;
				}

				return {
					callback : callback,
					ref : ref || null,
					priority : parseInt(priority) || 100,
					once : (once == true)
				};
			}

			//register a subscriber, handle if scope, return unsubscriber function
			function registerSubscriber (topic, subscriber) {
				if (subscriber && !angular.isEmpty(subscriber)) {
					if(!map[topic]) {
						map[topic] = [];
					}

					//if ref is scope, register destroy
					var ref = subscriber.ref;
					if (angular.isScope(ref)) {
						ref.$on('$destroy', function() {
							destroy(parseRefId(ref));
						});
					}

					map[topic].push(subscriber);

					return angular.once(function () {
						unregisterSubscriber(topic, subscriber);
					});
				} else {
					return angular.noop;
				}
			}

			function unregisterSubscriber(topic, subscriber) {
				var removed = angular.remove(map[topic], function (sub) {
					return angular.equals(sub, subscriber);
				});
				return removed.length > 0;
			}

			//special handling for scopes
			function parseRefId(ref) {
				return angular.isScope(ref) ? ref.$id : ref;
			}

			/*********
			 Functions
			 **********/

			var logListeners = function pubsub_logListeners () {
				Debugger.log('LISTENERS:');
				angular.forEach(map, function (val, key) {
					Debugger.log(key);
					Debugger.table(map[key])
				});
			};

			/**
			 @name PubSub.trigger
			 @description
			 Publish some data on a topic
			 @param topic {string} channel to publish on, can be multiple space-separated
			 @param args {*}  Array of arguments to apply to callback. If you want to pass a single arg, you don't need to wrap it in an array.
			 */
			var trigger = function pubsub_trigger (topic, args) {

				//ensure arguments are array
				//undefined is not a JSON value, so means nothing was sent over the socket
				if (angular.isUndefined(args)) {
					args = null;
				} else {
					args = [args];
				}

				//loop through each passed topic
				angular.forEach(splitTopics(topic), function (current) {
					//loop through each subscriber
					Debugger.log('Publish on ' + current, args);
					if (checkTopicHasSubs(current)) {

						angular.forEach( $filter('orderBy')(map[current], 'priority') , function (subscriber, index) {
							//future - avoid $safeApply
							$rootScope.$safeApply(function() {
								subscriber.callback.apply(subscriber.ref, args);
							});

							if (subscriber.once == true) {
								map[current].splice(index, 1);
							}
						});
					}
				});
			};

			/**
			 @name PubSub.reject
			 @description
			 Cancel a callback, publish null on the topic
			 @param topic {string} channel to publish on, can be multiple space-separated
			*/
			var reject = function pubsub_reject (topic) {
				angular.forEach(splitTopics(topic), function (current) {
					Debugger.log('Reject on ' + current);
					angular.forEach(map[current], function (subscriber, index) {
						//future - avoid $safeApply
						$rootScope.$safeApply(function() {
							subscriber.callback.apply(subscriber.ref, null);
						});

						if (subscriber.once == true) {
							map[current].splice(index, 1);
						}
					});
				});
			};

			/**
			 @name PubSub.on
			 @description
			 Register a callback to a topic

			 @param topic {string} channel to subscribe to. Can pass in multiple, space-separated.
			 @param callback {function} handler event. Will be called on publish event to given topic, passed args from publish
			 @param ref {string} context for this, or $scope (will automatically set up $destroy listener), or reference ID to be used in PubSub.destroy()
			 @param priority {number} Priority to run function at. default is 100. Lower gets priority.
			 @param one {boolean} Flag to run the callback only once
			 @return {Function} handle to pass into unsubscribe. If multiple events are passed in, an array is returned.
			 */
			//note - ref is also context for this in callback
			var on = function pubsub_on (topic, callback, ref, priority, one) {
				ref = ref || null;
				one = one == true;

				if (one) {
					callback = angular.once(callback);
				}

				var unsubscribers = [];

				angular.forEach(splitTopics(topic), function (current) {
					unsubscribers.push(
						registerSubscriber(current,
							createSubscriber(callback, ref, priority, one)
						)
					);
				});

				return (function pubsub_unsubscribe() {
					angular.forEach(unsubscribers, function (handle) {
						handle();
					});
				});
			};

			/**
			 * @name PubSub.once
			 * @description
			 * Register a callback to a topic only a single time
			 * @param topic {string}
			 * @param callback {function}
			 * @param ref {string} context for this, or $scope (will automatically set up $destroy listener), or reference ID to be used in PubSub.destroy()
			 * @param priority {number} Priority to run function at. default is 100. Lower gets priority.
			 *
			 */
			var once = function pubsub_once(topic, callback, ref, priority) {
				on(topic, callback, ref, priority, true);
			};

			/**
			 @name PubSub.off
			 @param topic {string}
			 @param callback {function}
			 @description
			 Disconnect specific callback function from topic. Note, however, you usually would do this by running the function returned by on()

			 e.g. var handle = PubSub.on( "someTopic", function(){} );
			 handle(); //unsubscribe

			 */
			var off = function pubsub_off (topic, callback) {

				var removed = angular.remove(map[topic], function (subscriber) {
					return angular.equals(subscriber.callback, callback);
				});

				return (removed.length > 0);
			};

			/**
			 * @name PubSub.destroy
			 *
			 * @param ref {string|object} Remove all listeners with a given reference (same object passed on on())
			 *
			 * @description
			 * Removes all listeners for an associated reference. When pass a $scope to on(), destroy will automatically be set up on $scope.$destroy, so this is unneeded.
			 */
			var destroy = function pubsub_destroy(ref) {
				angular.forEach(map, function (topicArray, key) {
					angular.remove(topicArray, function (subscriber) {
						return angular.equals(parseRefId(ref), parseRefId(subscriber.ref));
					});
				});
			};

			/**
			 * @name PubSub.clear
			 *
			 * @param topic {string} Topic to clear, can be space-separated list
			 *
			 * @description
			 * Clears all listeners in a topic.
			 * Clears all listeners if "all" is passed.
			 * Does nothing if no topic passed.
			 */
			//remove all subscribers for given topic(s)
			var clear = function pubsub_clear (topic) {
				angular.forEach(splitTopics(topic), function (val, key) {
					map[key].length = 0;
				});
			};

			return {
				logListeners: logListeners,
				trigger: trigger,
				on: on,
				once: once,
				off: off,
				destroy: destroy,
				reject : reject,
				clear: clear
			}
		}
	});
