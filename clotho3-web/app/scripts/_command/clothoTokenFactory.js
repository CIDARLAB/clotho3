angular.module('clotho.tokenizer')
/**
 * @name clothoTokenFactory
 *
 * @description
 * ClothoTokens are essentially wrappers for clotho sharables, or strings.
 * They expect the fields minimally of name, id, schema to be a sharable, or just a string for other keywords
 */
	.factory('clothoTokenFactory', function (Clotho) {

		//pass UUID to make object, or just pass value as string
		function ClothoToken (sharable) {
			var self = this;
			self.model = sharable;

			if (this.isSharable()) {
				self.fullSharablePromise = Clotho.get(self.model.id, {mute : true}).then(function (data) {
					self.fullSharable = data;
				});
			}
		}

		ClothoToken.prototype.readable = function () {
			return this.model.name || this.model;
		};

		ClothoToken.prototype.id = function () {
			return this.model.id || null;
		};

		ClothoToken.prototype.isAmbiguous = function () {
			return angular.isArray(this.model);
		};

		ClothoToken.prototype.isSharable = function () {
			return !this.isAmbiguous() && angular.isDefined(this.model.id);
		};

		return ClothoToken;

	});