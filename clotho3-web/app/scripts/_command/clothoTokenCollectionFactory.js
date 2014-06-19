angular.module('clotho.tokenizer')
/**
 * @name clothoTokenCollectionFactory
 *
 * @description
 * Object to handle a collection of ClothoTokens
 */
	.factory('clothoTokenCollectionFactory', function (clothoTokenFactory) {

		function ClothoTokenCollection (startingTokens) {

			var self = this;
			this.tokens = [];
			this.currentTokenIndex = -1;

			if (angular.isArray(startingTokens) || angular.isObject(startingTokens)) {
				angular.forEach(startingTokens, function (token) {
					self.addToken(token);
				});
			}
		}

		//add a token, pass arguments through to clothoTokenFactory
		ClothoTokenCollection.prototype.addToken = function (sharable) {
			this.tokens.push(new clothoTokenFactory(sharable));
		};

		ClothoTokenCollection.prototype.inRange = function (index) {
			return index > -1 && index < this.tokens.length;
		};

		//get token at given index
		ClothoTokenCollection.prototype.getToken = function (index) {
			return this.tokens[index];
		};

		//return index of token
		ClothoTokenCollection.prototype.indexOf = function (token) {
			return this.tokens.indexOf(token);
		};

		//remove token at given index, return it if removed, otherwise false
		ClothoTokenCollection.prototype.removeToken = function (index) {
			if (this.inRange(index)) {
				return this.tokens.splice(index, 1);
			} else {
				return false;
			}
		};

		//remove all tokens
		ClothoTokenCollection.prototype.removeAll = function () {
			this.tokens.length = 0;
		};

		//remove active token if set and return it, otherwise return false
		ClothoTokenCollection.prototype.removeActiveToken = function () {
			if (this.isActive()) {
				var toReturn =  this.removeToken(this.currentTokenIndex);
				this.unsetActive();
				return toReturn;
			} else {
				return false;
			}
		};

		//set token at given index to be active
		ClothoTokenCollection.prototype.setActive = function (index) {
			if (this.inRange(index)) {
				this.currentTokenIndex = index;
				return index;
			} else {
				return false;
			}
		};

		ClothoTokenCollection.prototype.toggleActive = function (index) {
			this[this.isActive(index) ? 'unsetActive' : 'setActive'](index);
		};

		//set token at last position to be active
		ClothoTokenCollection.prototype.setLastActive = function (index) {
			this.setActive(this.tokens.length - 1);
		};

		//set previous token active, based on current active, otherwise last
		ClothoTokenCollection.prototype.setPrevActive = function (index) {
			this.currentTokenIndex = (this.currentTokenIndex > 0 ? this.currentTokenIndex : this.tokens.length) - 1;
		};

		//set next token active, based on current active, otherwise first
		ClothoTokenCollection.prototype.setNextActive = function (index) {
			this.currentTokenIndex = (this.currentTokenIndex + 1) % this.tokens.length;
		};

		//unset active token
		ClothoTokenCollection.prototype.unsetActive = function (index) {
			this.currentTokenIndex = -1;
		};

		//check if token is active, at index when passed, otherwise if any is active
		ClothoTokenCollection.prototype.isActive = function (index) {
			if (angular.isDefined(index)) {
				return this.currentTokenIndex == index;
			} else {
				return this.currentTokenIndex > -1;
			}
		};

		ClothoTokenCollection.prototype.whichActive = function () {
			return this.currentTokenIndex;
		};

		ClothoTokenCollection.prototype.isLastActive = function () {
			return (this.currentTokenIndex === (this.tokens.length - 1))
		};

		return ClothoTokenCollection;
	});